/// VCF file parsing using noodles-vcf.
///
/// Handles single-sample and multi-sample VCFs, multi-allelic decomposition,
/// and extraction of INFO/FORMAT fields as JSON.
use std::io::BufRead;
use std::path::Path;

use noodles::bgzf;
use noodles::vcf;
use noodles::vcf::variant::record::samples::series::value::genotype::Phasing;
use serde_json::{json, Map, Value};

use crate::assembly::{Assembly, ChrIndex};
use crate::uvid128::Uvid128;

/// A site-level record (shared across samples).
#[derive(Debug, Clone)]
pub struct SiteRecord {
    pub uvid: Uvid128,
    pub qual: Option<f32>,
    pub filter: String,
    pub info: Value,
    pub multiallelic: bool,
}

/// A sample-level record (per-sample genotype data).
#[derive(Debug, Clone)]
pub struct SampleRecord {
    pub uvid: Uvid128,
    pub allele1: Option<u8>,
    pub allele2: Option<u8>,
    pub phased: bool,
    pub dp: Option<u16>,
    pub gq: Option<u8>,
    pub format_extra: Value,
}

/// Result of parsing a VCF file: site records and per-sample records.
#[derive(Debug)]
pub struct ParsedVcf {
    pub sample_names: Vec<String>,
    pub sites: Vec<SiteRecord>,
    /// Per-sample records, indexed by sample index (matching sample_names order).
    pub samples: Vec<Vec<SampleRecord>>,
    pub raw_header: String,
    pub info_definitions: Value,
    pub format_definitions: Value,
    pub contigs: Value,
}

/// Parse a VCF file (plain or bgzf-compressed) and decompose multi-allelic sites.
pub fn parse_vcf(path: &Path, assembly: Assembly) -> Result<ParsedVcf, Box<dyn std::error::Error>> {
    let is_bgzf = {
        let mut peek = std::fs::File::open(path)?;
        let mut magic = [0u8; 2];
        use std::io::Read;
        match peek.read_exact(&mut magic) {
            Ok(()) => magic == [0x1f, 0x8b],
            Err(_) => false,
        }
    };

    if is_bgzf {
        let file = std::fs::File::open(path)?;
        let decoder = bgzf::io::Reader::new(file);
        let buf_reader = std::io::BufReader::new(decoder);
        parse_vcf_from_reader(buf_reader, assembly)
    } else {
        let file = std::fs::File::open(path)?;
        let buf_reader = std::io::BufReader::new(file);
        parse_vcf_from_reader(buf_reader, assembly)
    }
}

fn parse_vcf_from_reader<R: BufRead>(
    reader: R,
    assembly: Assembly,
) -> Result<ParsedVcf, Box<dyn std::error::Error>> {
    // Use the noodles VCF reader to handle header + record parsing properly
    let mut vcf_reader = vcf::io::Reader::new(reader);

    // Read header
    let header = vcf_reader.read_header()?;
    // Serialize header back to string via the writer
    let raw_header = {
        let mut buf = Vec::new();
        let mut writer = vcf::io::Writer::new(&mut buf);
        writer.write_header(&header)?;
        String::from_utf8_lossy(&buf).to_string()
    };

    // Extract sample names from header
    let sample_names: Vec<String> = header
        .sample_names()
        .iter()
        .map(|s| s.to_string())
        .collect();

    // Extract metadata
    let info_definitions = extract_info_definitions(&header);
    let format_definitions = extract_format_definitions(&header);
    let contigs = extract_contig_definitions(&header);

    let num_samples = sample_names.len();
    let mut sites: Vec<SiteRecord> = Vec::new();
    let mut samples: Vec<Vec<SampleRecord>> = (0..num_samples).map(|_| Vec::new()).collect();

    // Parse records using RecordBuf
    let mut record_buf = vcf::variant::RecordBuf::default();
    loop {
        match vcf_reader.read_record_buf(&header, &mut record_buf) {
            Ok(0) => break, // EOF
            Ok(_) => {}
            Err(e) => return Err(e.into()),
        }

        let chr_name = record_buf.reference_sequence_name().to_string();
        let chr = match ChrIndex::from_name(&chr_name) {
            Some(c) => c,
            None => continue,
        };

        let pos: usize = record_buf.variant_start().map(|p| p.get()).unwrap_or(0);
        if pos == 0 || pos > u32::MAX as usize {
            continue;
        }
        let pos = pos as u32;

        let ref_seq = record_buf.reference_bases().to_string();

        // Collect ALT alleles via AsRef on the concrete AlternateBases type
        let alt_alleles: Vec<String> = record_buf
            .alternate_bases()
            .as_ref()
            .iter()
            .map(|a| a.to_string())
            .collect();

        let is_multiallelic = alt_alleles.len() > 1;

        // quality_score() on RecordBuf returns Option<f32> directly
        let qual = record_buf.quality_score();

        let filter_str = format_filter(&record_buf);
        let info_json = format_info(&record_buf);

        // For each ALT allele, create a site record and sample records
        let alts_to_process = if alt_alleles.is_empty() {
            vec![".".to_string()]
        } else {
            alt_alleles.clone()
        };

        for (alt_idx, alt_allele) in alts_to_process.iter().enumerate() {
            let alt_bytes = if alt_allele == "." || alt_allele == "*" {
                b".".to_vec()
            } else {
                alt_allele.as_bytes().to_vec()
            };

            let uvid = match Uvid128::encode(chr, pos, ref_seq.as_bytes(), &alt_bytes, assembly) {
                Some(u) => u,
                None => continue,
            };

            sites.push(SiteRecord {
                uvid,
                qual,
                filter: filter_str.clone(),
                info: info_json.clone(),
                multiallelic: is_multiallelic,
            });

            // Process each sample's genotype at this site
            for (sample_idx, _sample_name) in sample_names.iter().enumerate() {
                let (allele1, allele2, phased, dp, gq, format_extra) =
                    extract_sample_data(&record_buf, sample_idx, alt_idx);

                // Skip homref (both alleles are 0)
                let is_homref = allele1 == Some(0) && allele2 == Some(0);

                if is_homref {
                    continue;
                }

                samples[sample_idx].push(SampleRecord {
                    uvid,
                    allele1,
                    allele2,
                    phased,
                    dp,
                    gq,
                    format_extra,
                });
            }
        }
    }

    Ok(ParsedVcf {
        sample_names,
        sites,
        samples,
        raw_header,
        info_definitions,
        format_definitions,
        contigs,
    })
}

/// Extract all relevant data from a sample at a given site.
fn extract_sample_data(
    record: &vcf::variant::RecordBuf,
    sample_idx: usize,
    _alt_idx: usize,
) -> (Option<u8>, Option<u8>, bool, Option<u16>, Option<u8>, Value) {
    use vcf::variant::record_buf::samples::sample::Value as SampleValue;

    // Default no-call result
    let mut allele1: Option<u8> = None;
    let mut allele2: Option<u8> = None;
    let mut phased = false;
    let mut dp: Option<u16> = None;
    let mut gq: Option<u8> = None;
    let mut extra = Map::new();

    // Get the sample by index
    let sample = match record.samples().get_index(sample_idx) {
        Some(s) => s,
        None => return (None, None, false, None, None, Value::Null),
    };

    // Iterate over keys and values together
    let keys = record.samples().keys();
    let values = sample.values();

    for (key_idx, key) in keys.as_ref().iter().enumerate() {
        let value = match values.get(key_idx) {
            Some(Some(v)) => v,
            _ => continue,
        };

        match key.as_str() {
            "GT" => {
                // Parse genotype
                if let SampleValue::Genotype(gt) = value {
                    let alleles = gt.as_ref(); // &[Allele]
                    if let Some(first) = alleles.first() {
                        allele1 = first.position().map(|p| p as u8);
                    }
                    if alleles.len() > 1 {
                        if let Some(second) = alleles.get(1) {
                            allele2 = second.position().map(|p| p as u8);
                            phased = second.phasing() == Phasing::Phased;
                        }
                    }
                }
            }
            "DP" => {
                if let SampleValue::Integer(v) = value {
                    dp = Some(*v as u16);
                }
            }
            "GQ" => {
                if let SampleValue::Integer(v) = value {
                    gq = Some(*v as u8);
                }
            }
            _ => {
                // Store other FORMAT fields in format_extra
                extra.insert(key.clone(), sample_value_to_json(value));
            }
        }
    }

    let format_extra = if extra.is_empty() {
        Value::Null
    } else {
        Value::Object(extra)
    };

    (allele1, allele2, phased, dp, gq, format_extra)
}

/// Convert a sample field value to JSON.
fn sample_value_to_json(value: &vcf::variant::record_buf::samples::sample::Value) -> Value {
    use vcf::variant::record_buf::samples::sample::Value as V;
    match value {
        V::Integer(n) => json!(n),
        V::Float(f) => json!(f),
        V::Character(c) => json!(c.to_string()),
        V::String(s) => json!(s),
        V::Genotype(gt) => json!(format!("{:?}", gt)),
        V::Array(arr) => {
            use vcf::variant::record_buf::samples::sample::value::Array as A;
            match arr {
                A::Integer(vals) => Value::Array(
                    vals.iter()
                        .map(|v| match v {
                            Some(n) => json!(n),
                            None => Value::Null,
                        })
                        .collect(),
                ),
                A::Float(vals) => Value::Array(
                    vals.iter()
                        .map(|v| match v {
                            Some(f) => json!(f),
                            None => Value::Null,
                        })
                        .collect(),
                ),
                A::Character(vals) => Value::Array(
                    vals.iter()
                        .map(|v| match v {
                            Some(c) => json!(c.to_string()),
                            None => Value::Null,
                        })
                        .collect(),
                ),
                A::String(vals) => Value::Array(
                    vals.iter()
                        .map(|v| match v {
                            Some(s) => json!(s),
                            None => Value::Null,
                        })
                        .collect(),
                ),
            }
        }
    }
}

/// Extract INFO field definitions from the VCF header as JSON.
fn extract_info_definitions(header: &vcf::Header) -> Value {
    let defs: Vec<Value> = header
        .infos()
        .iter()
        .map(|(id, info)| {
            json!({
                "id": id.to_string(),
                "number": format!("{:?}", info.number()),
                "type": format!("{:?}", info.ty()),
                "description": info.description(),
            })
        })
        .collect();
    Value::Array(defs)
}

/// Extract FORMAT field definitions from the VCF header as JSON.
fn extract_format_definitions(header: &vcf::Header) -> Value {
    let defs: Vec<Value> = header
        .formats()
        .iter()
        .map(|(id, format)| {
            json!({
                "id": id.to_string(),
                "number": format!("{:?}", format.number()),
                "type": format!("{:?}", format.ty()),
                "description": format.description(),
            })
        })
        .collect();
    Value::Array(defs)
}

/// Extract contig definitions from the VCF header as JSON.
fn extract_contig_definitions(header: &vcf::Header) -> Value {
    let defs: Vec<Value> = header
        .contigs()
        .iter()
        .map(|(id, contig)| {
            let mut map = Map::new();
            map.insert("id".to_string(), json!(id.to_string()));
            if let Some(length) = contig.length() {
                map.insert("length".to_string(), json!(length));
            }
            Value::Object(map)
        })
        .collect();
    Value::Array(defs)
}

/// Format the FILTER column as a string.
fn format_filter(record: &vcf::variant::RecordBuf) -> String {
    let filters = record.filters();
    if filters.as_ref().is_empty() {
        ".".to_string()
    } else {
        let filter_strs: Vec<&str> = filters.as_ref().iter().map(|f| f.as_str()).collect();
        filter_strs.join(";")
    }
}

/// Format INFO fields as a JSON object.
fn format_info(record: &vcf::variant::RecordBuf) -> Value {
    let info = record.info();
    let mut map = Map::new();

    for (key, value) in info.as_ref() {
        let key_str = key.to_string();
        match value {
            Some(v) => {
                map.insert(key_str, info_value_to_json(v));
            }
            None => {
                map.insert(key_str, Value::Bool(true));
            }
        }
    }

    Value::Object(map)
}

/// Convert a noodles INFO field value to JSON.
fn info_value_to_json(value: &vcf::variant::record_buf::info::field::Value) -> Value {
    use vcf::variant::record_buf::info::field::Value as V;
    match value {
        V::Integer(n) => json!(n),
        V::Float(f) => json!(f),
        V::Character(c) => json!(c.to_string()),
        V::String(s) => json!(s),
        V::Flag => Value::Bool(true),
        V::Array(arr) => {
            use vcf::variant::record_buf::info::field::value::Array as A;
            match arr {
                A::Integer(vals) => Value::Array(
                    vals.iter()
                        .map(|v| match v {
                            Some(n) => json!(n),
                            None => Value::Null,
                        })
                        .collect(),
                ),
                A::Float(vals) => Value::Array(
                    vals.iter()
                        .map(|v| match v {
                            Some(f) => json!(f),
                            None => Value::Null,
                        })
                        .collect(),
                ),
                A::Character(vals) => Value::Array(
                    vals.iter()
                        .map(|v| match v {
                            Some(c) => json!(c.to_string()),
                            None => Value::Null,
                        })
                        .collect(),
                ),
                A::String(vals) => Value::Array(
                    vals.iter()
                        .map(|v| match v {
                            Some(s) => json!(s),
                            None => Value::Null,
                        })
                        .collect(),
                ),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_info_value_to_json_integer() {
        use vcf::variant::record_buf::info::field::Value as V;
        let val = V::Integer(42);
        assert_eq!(info_value_to_json(&val), json!(42));
    }

    #[test]
    fn test_info_value_to_json_float() {
        use vcf::variant::record_buf::info::field::Value as V;
        let val = V::Float(0.5);
        assert_eq!(info_value_to_json(&val), json!(0.5));
    }

    #[test]
    fn test_info_value_to_json_string() {
        use vcf::variant::record_buf::info::field::Value as V;
        let val = V::String("hello".to_string());
        assert_eq!(info_value_to_json(&val), json!("hello"));
    }

    #[test]
    fn test_info_value_to_json_flag() {
        use vcf::variant::record_buf::info::field::Value as V;
        let val = V::Flag;
        assert_eq!(info_value_to_json(&val), json!(true));
    }
}
