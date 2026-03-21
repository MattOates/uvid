/// DuckDB-backed storage for .uvid collection files.
///
/// A `.uvid` file is a DuckDB database containing:
/// - `<file>__sites` tables: variant site data shared across samples
/// - `<file>__<sample>` tables: per-sample genotype data
/// - `_meta` table: provenance and VCF header metadata
use std::path::Path;

use duckdb::{params, Connection};
use serde_json::Value;

use crate::uvid128::Uvid128;
use crate::vcf::ParsedVcf;

/// A handle to a .uvid collection file (DuckDB database).
pub struct UvidStore {
    conn: Connection,
}

impl UvidStore {
    /// Open or create a .uvid collection file.
    pub fn open(path: &Path) -> Result<Self, Box<dyn std::error::Error>> {
        let conn = Connection::open(path)?;
        let store = UvidStore { conn };
        store.ensure_meta_table()?;
        // Install and load JSON extension
        store.conn.execute_batch("INSTALL json; LOAD json;").ok(); // Ignore if already installed
        Ok(store)
    }

    /// Open an in-memory store (for testing).
    pub fn open_in_memory() -> Result<Self, Box<dyn std::error::Error>> {
        let conn = Connection::open_in_memory()?;
        let store = UvidStore { conn };
        store.ensure_meta_table()?;
        store.conn.execute_batch("INSTALL json; LOAD json;").ok();
        Ok(store)
    }

    /// Create the _meta table if it doesn't exist.
    fn ensure_meta_table(&self) -> Result<(), Box<dyn std::error::Error>> {
        self.conn.execute_batch(
            "CREATE TABLE IF NOT EXISTS _meta (
                table_name VARCHAR PRIMARY KEY,
                table_type VARCHAR NOT NULL,
                source_file VARCHAR NOT NULL,
                sample_name VARCHAR,
                vcf_header TEXT NOT NULL,
                info_definitions JSON,
                format_definitions JSON,
                contigs JSON,
                added_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )",
        )?;
        Ok(())
    }

    /// Derive the table name prefix from a VCF filename.
    ///
    /// Strips the extension(s) from the filename:
    /// - `sample.vcf` -> `sample`
    /// - `sample.vcf.gz` -> `sample`
    /// - `/path/to/sample.vcf.gz` -> `sample`
    pub fn table_prefix(vcf_path: &Path) -> String {
        let stem = vcf_path
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();

        // Strip .vcf.gz, .vcf, .bcf, .gz
        let stripped = stem
            .strip_suffix(".vcf.gz")
            .or_else(|| stem.strip_suffix(".vcf"))
            .or_else(|| stem.strip_suffix(".bcf"))
            .or_else(|| stem.strip_suffix(".gz"))
            .unwrap_or(&stem);

        stripped.to_string()
    }

    /// Add a parsed VCF to the collection.
    ///
    /// Creates `<prefix>__sites` and `<prefix>__<sample>` tables.
    pub fn add_vcf(
        &self,
        prefix: &str,
        parsed: &ParsedVcf,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let sites_table = format!("{}__sites", prefix);

        // Create sites table
        self.conn.execute_batch(&format!(
            "CREATE TABLE IF NOT EXISTS \"{}\" (
                uvid HUGEINT PRIMARY KEY,
                qual FLOAT,
                filter VARCHAR,
                info JSON,
                multiallelic BOOLEAN NOT NULL DEFAULT FALSE
            )",
            sites_table,
        ))?;

        // Bulk insert sites
        {
            let mut appender = self.conn.appender(&sites_table)?;
            for site in &parsed.sites {
                let uvid_i128 = site.uvid.as_u128() as i128;
                let info_str = site.info.to_string();
                appender.append_row(params![
                    uvid_i128,
                    site.qual,
                    site.filter,
                    info_str,
                    site.multiallelic,
                ])?;
            }
            appender.flush()?;
        }

        // Register sites table in _meta
        self.conn.execute(
            "INSERT OR REPLACE INTO _meta (table_name, table_type, source_file, sample_name, vcf_header, info_definitions, format_definitions, contigs)
             VALUES (?, 'sites', ?, NULL, ?, ?, ?, ?)",
            params![
                sites_table,
                prefix,
                parsed.raw_header,
                parsed.info_definitions.to_string(),
                parsed.format_definitions.to_string(),
                parsed.contigs.to_string(),
            ],
        )?;

        // Create and populate per-sample tables
        for (sample_idx, sample_name) in parsed.sample_names.iter().enumerate() {
            let sample_table = format!("{}__{}", prefix, sample_name);

            self.conn.execute_batch(&format!(
                "CREATE TABLE IF NOT EXISTS \"{}\" (
                    uvid HUGEINT PRIMARY KEY,
                    allele1 UTINYINT,
                    allele2 UTINYINT,
                    phased BOOLEAN,
                    dp USMALLINT,
                    gq UTINYINT,
                    format_extra JSON
                )",
                sample_table,
            ))?;

            // Bulk insert sample records
            {
                let mut appender = self.conn.appender(&sample_table)?;
                for record in &parsed.samples[sample_idx] {
                    let uvid_i128 = record.uvid.as_u128() as i128;
                    let format_extra_str = if record.format_extra.is_null() {
                        None
                    } else {
                        Some(record.format_extra.to_string())
                    };
                    appender.append_row(params![
                        uvid_i128,
                        record.allele1,
                        record.allele2,
                        record.phased,
                        record.dp,
                        record.gq,
                        format_extra_str,
                    ])?;
                }
                appender.flush()?;
            }

            // Register sample table in _meta
            self.conn.execute(
                "INSERT OR REPLACE INTO _meta (table_name, table_type, source_file, sample_name, vcf_header, info_definitions, format_definitions, contigs)
                 VALUES (?, 'sample', ?, ?, ?, ?, ?, ?)",
                params![
                    sample_table,
                    prefix,
                    sample_name,
                    parsed.raw_header,
                    parsed.info_definitions.to_string(),
                    parsed.format_definitions.to_string(),
                    parsed.contigs.to_string(),
                ],
            )?;
        }

        Ok(())
    }

    /// List all sample table names in the collection.
    pub fn list_samples(
        &self,
    ) -> Result<Vec<(String, String, String)>, Box<dyn std::error::Error>> {
        let mut stmt = self.conn.prepare(
            "SELECT table_name, source_file, sample_name FROM _meta WHERE table_type = 'sample' ORDER BY table_name"
        )?;

        let rows = stmt.query_map([], |row| {
            Ok((
                row.get::<_, String>(0)?,
                row.get::<_, String>(1)?,
                row.get::<_, String>(2)?,
            ))
        })?;

        let mut result = Vec::new();
        for row in rows {
            result.push(row?);
        }
        Ok(result)
    }

    /// List all source files in the collection.
    pub fn list_sources(&self) -> Result<Vec<String>, Box<dyn std::error::Error>> {
        let mut stmt = self
            .conn
            .prepare("SELECT DISTINCT source_file FROM _meta ORDER BY source_file")?;

        let rows = stmt.query_map([], |row| row.get::<_, String>(0))?;

        let mut result = Vec::new();
        for row in rows {
            result.push(row?);
        }
        Ok(result)
    }

    /// Query variants in a genomic region from a sample table.
    pub fn search_region(
        &self,
        table_name: &str,
        chr: crate::assembly::ChrIndex,
        start_pos: u32,
        end_pos: u32,
    ) -> Result<Vec<SearchResult>, Box<dyn std::error::Error>> {
        let (lower, upper) = Uvid128::range(chr, start_pos, end_pos);

        let query = format!(
            "SELECT s.uvid, s.allele1, s.allele2, s.phased, s.dp, s.gq, s.format_extra,
                    site.qual, site.filter, site.info, site.multiallelic
             FROM \"{}\" s
             LEFT JOIN \"{}\" site ON s.uvid = site.uvid
             WHERE s.uvid BETWEEN ? AND ?
             ORDER BY s.uvid",
            table_name,
            table_name.replace(
                &format!("__{}", table_name.rsplit("__").next().unwrap_or("")),
                "__sites",
            ),
        );

        // This is a simplified approach; the join table name derivation
        // should be done more robustly in practice
        let lower_i128 = lower.as_u128() as i128;
        let upper_i128 = upper.as_u128() as i128;

        let mut stmt = self.conn.prepare(&query)?;
        let rows = stmt.query_map(params![lower_i128, upper_i128], |row| {
            let uvid_i128: i128 = row.get(0)?;
            Ok(SearchResult {
                uvid: Uvid128::from_u128(uvid_i128 as u128),
                allele1: row.get(1)?,
                allele2: row.get(2)?,
                phased: row.get(3)?,
                dp: row.get(4)?,
                gq: row.get(5)?,
                format_extra: row
                    .get::<_, Option<String>>(6)?
                    .map(|s| serde_json::from_str(&s).unwrap_or(Value::Null))
                    .unwrap_or(Value::Null),
                qual: row.get(7)?,
                filter: row.get(8)?,
                info: row
                    .get::<_, Option<String>>(9)?
                    .map(|s| serde_json::from_str(&s).unwrap_or(Value::Null))
                    .unwrap_or(Value::Null),
                multiallelic: row.get(10)?,
            })
        })?;

        let mut results = Vec::new();
        for row in rows {
            results.push(row?);
        }
        Ok(results)
    }

    /// Execute a raw SQL query (for advanced use / debugging).
    pub fn execute_sql(&self, sql: &str) -> Result<(), Box<dyn std::error::Error>> {
        self.conn.execute_batch(sql)?;
        Ok(())
    }

    /// Get a reference to the underlying DuckDB connection.
    pub fn connection(&self) -> &Connection {
        &self.conn
    }
}

/// Result of a search query combining site and sample data.
#[derive(Debug, Clone)]
pub struct SearchResult {
    pub uvid: Uvid128,
    pub allele1: Option<u8>,
    pub allele2: Option<u8>,
    pub phased: Option<bool>,
    pub dp: Option<u16>,
    pub gq: Option<u8>,
    pub format_extra: Value,
    pub qual: Option<f32>,
    pub filter: Option<String>,
    pub info: Value,
    pub multiallelic: Option<bool>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_table_prefix() {
        assert_eq!(UvidStore::table_prefix(Path::new("sample.vcf")), "sample");
        assert_eq!(
            UvidStore::table_prefix(Path::new("sample.vcf.gz")),
            "sample"
        );
        assert_eq!(
            UvidStore::table_prefix(Path::new("/path/to/sample.vcf.gz")),
            "sample"
        );
        assert_eq!(UvidStore::table_prefix(Path::new("my_data.bcf")), "my_data");
        assert_eq!(UvidStore::table_prefix(Path::new("noext")), "noext");
    }

    #[test]
    fn test_open_in_memory() {
        let store = UvidStore::open_in_memory().unwrap();
        let samples = store.list_samples().unwrap();
        assert!(samples.is_empty());
    }

    #[test]
    fn test_meta_table_creation() {
        let store = UvidStore::open_in_memory().unwrap();
        let sources = store.list_sources().unwrap();
        assert!(sources.is_empty());
    }
}
