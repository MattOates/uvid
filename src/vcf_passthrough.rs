/// VCF passthrough: read a VCF record-by-record, replace the ID column with
/// UVID identifiers (or optionally UUIDv5 representations), and write the
/// modified VCF to stdout or a file.
///
/// Uses line-based processing for maximum throughput — only the ID column
/// (column 2, 0-indexed) is modified; every other byte passes through unchanged.
use std::fmt;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use noodles::bgzf;

use crate::assembly::{Assembly, ChrIndex};
use crate::uvid128::Uvid128;

// ---------------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------------

/// Errors that can occur during VCF passthrough processing.
#[derive(Debug)]
pub enum VcfPassthroughError {
    /// Assembly could not be detected from the VCF header and no override was provided.
    AssemblyNotDetected,
    /// An I/O error occurred while reading or writing.
    Io(io::Error),
}

impl fmt::Display for VcfPassthroughError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VcfPassthroughError::AssemblyNotDetected => write!(
                f,
                "Assembly could not be detected from VCF header. \
                 Use --assembly to specify GRCh37 or GRCh38."
            ),
            VcfPassthroughError::Io(e) => write!(f, "I/O error: {}", e),
        }
    }
}

impl std::error::Error for VcfPassthroughError {}

impl From<io::Error> for VcfPassthroughError {
    fn from(e: io::Error) -> Self {
        VcfPassthroughError::Io(e)
    }
}

// ---------------------------------------------------------------------------
// Assembly detection
// ---------------------------------------------------------------------------

/// Try to detect the genome assembly from VCF header lines.
///
/// Scans for patterns in:
/// - `##reference=...`  (e.g. `##reference=GRCh38` or `##reference=file:///path/to/hg19.fa`)
/// - `##assembly=...`   (e.g. `##assembly=GRCh38`)
/// - `##contig=<...,assembly=...>` (e.g. `##contig=<ID=chr1,length=248956422,assembly=GRCh38>`)
///
/// Returns `Some(Assembly)` on the first unambiguous match, `None` otherwise.
pub fn detect_assembly_from_header(header_lines: &[String]) -> Option<Assembly> {
    for line in header_lines {
        if let Some(rest) = line.strip_prefix("##reference=") {
            if let Some(asm) = match_assembly_pattern(rest) {
                return Some(asm);
            }
        } else if let Some(rest) = line.strip_prefix("##assembly=") {
            if let Some(asm) = match_assembly_pattern(rest) {
                return Some(asm);
            }
        } else if line.starts_with("##contig=<") {
            // Look for assembly= key inside the angle-bracket body
            if let Some(start) = line.find("assembly=") {
                let rest = &line[start + "assembly=".len()..];
                // Value ends at ',' or '>'
                let end = rest.find([',', '>']).unwrap_or(rest.len());
                let value = &rest[..end];
                if let Some(asm) = match_assembly_pattern(value) {
                    return Some(asm);
                }
            }
        }
    }
    None
}

/// Check if a string contains a recognisable assembly identifier.
///
/// Matches case-insensitively anywhere within the string:
/// `GRCh37`, `GRCh38`, `hg19`, `hg38`.
fn match_assembly_pattern(s: &str) -> Option<Assembly> {
    let lower = s.to_lowercase();
    if lower.contains("grch38") || lower.contains("hg38") {
        Some(Assembly::GRCh38)
    } else if lower.contains("grch37") || lower.contains("hg19") {
        Some(Assembly::GRCh37)
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Core passthrough
// ---------------------------------------------------------------------------

/// Process a VCF file, replacing the ID column with UVID hex strings
/// (or UUIDv5 if `use_uuid` is true).
///
/// - `input`            — path to the input VCF (.vcf or .vcf.gz)
/// - `output`           — output path; `None` → write plain VCF to stdout.
///                         If the path ends in `.vcf.gz`, output is bgzf-compressed.
/// - `use_uuid`         — when true, emit UUIDv5 representation instead of UVID hex
/// - `assembly_override`— when `Some`, skip header detection and use this assembly
///
/// Returns the number of data records processed.
pub fn vcf_passthrough(
    input: &Path,
    output: Option<&Path>,
    use_uuid: bool,
    assembly_override: Option<Assembly>,
) -> Result<u64, VcfPassthroughError> {
    // Detect bgzf by magic bytes
    let is_bgzf = {
        let mut file = std::fs::File::open(input)?;
        let mut magic = [0u8; 2];
        use std::io::Read;
        match file.read_exact(&mut magic) {
            Ok(()) => magic == [0x1f, 0x8b],
            Err(_) => false,
        }
    };

    if is_bgzf {
        let file = std::fs::File::open(input)?;
        let decoder = bgzf::io::Reader::new(file);
        let reader = BufReader::new(decoder);
        passthrough_from_reader(reader, output, use_uuid, assembly_override)
    } else {
        let file = std::fs::File::open(input)?;
        let reader = BufReader::new(file);
        passthrough_from_reader(reader, output, use_uuid, assembly_override)
    }
}

/// Inner passthrough implementation generic over reader type.
fn passthrough_from_reader<R: BufRead>(
    mut reader: R,
    output: Option<&Path>,
    use_uuid: bool,
    assembly_override: Option<Assembly>,
) -> Result<u64, VcfPassthroughError> {
    // ---- Phase 1: read and pass through header lines, detect assembly ----
    let mut header_lines: Vec<String> = Vec::new();
    let mut first_data_line: Option<String> = None;
    let mut line_buf = String::new();

    loop {
        line_buf.clear();
        let n = reader.read_line(&mut line_buf)?;
        if n == 0 {
            break; // EOF before any data
        }
        let trimmed = line_buf.trim_end_matches(['\n', '\r']).to_string();
        if trimmed.starts_with('#') {
            header_lines.push(trimmed);
        } else {
            first_data_line = Some(trimmed);
            break;
        }
    }

    // Resolve assembly
    let assembly = match assembly_override {
        Some(asm) => asm,
        None => detect_assembly_from_header(&header_lines)
            .ok_or(VcfPassthroughError::AssemblyNotDetected)?,
    };

    // ---- Phase 2: open writer ----
    let wants_bgzf = output
        .map(|p| {
            p.to_string_lossy().ends_with(".vcf.gz") || p.to_string_lossy().ends_with(".vcf.bgz")
        })
        .unwrap_or(false);

    // We need to handle three cases: bgzf file, plain file, stdout.
    // Use an enum-dispatch approach to avoid Box<dyn Write>.
    enum OutputWriter {
        Bgzf(BufWriter<bgzf::io::Writer<std::fs::File>>),
        File(BufWriter<std::fs::File>),
        Stdout(BufWriter<io::Stdout>),
    }

    impl Write for OutputWriter {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            match self {
                OutputWriter::Bgzf(w) => w.write(buf),
                OutputWriter::File(w) => w.write(buf),
                OutputWriter::Stdout(w) => w.write(buf),
            }
        }
        fn flush(&mut self) -> io::Result<()> {
            match self {
                OutputWriter::Bgzf(w) => w.flush(),
                OutputWriter::File(w) => w.flush(),
                OutputWriter::Stdout(w) => w.flush(),
            }
        }
    }

    let mut writer: OutputWriter = if let Some(out_path) = output {
        if wants_bgzf {
            let file = std::fs::File::create(out_path)?;
            let bgzf_writer = bgzf::io::Writer::new(file);
            OutputWriter::Bgzf(BufWriter::new(bgzf_writer))
        } else {
            let file = std::fs::File::create(out_path)?;
            OutputWriter::File(BufWriter::new(file))
        }
    } else {
        OutputWriter::Stdout(BufWriter::new(io::stdout()))
    };

    // Write header lines
    for hline in &header_lines {
        writer.write_all(hline.as_bytes())?;
        writer.write_all(b"\n")?;
    }

    // ---- Phase 3: process data lines ----
    let mut count: u64 = 0;

    // Helper closure to process one data line
    let process_line = |line: &str, w: &mut OutputWriter| -> Result<(), VcfPassthroughError> {
        // Split into tab-separated fields. We need at least 5 columns
        // (CHROM, POS, ID, REF, ALT) to do anything useful.
        let fields: Vec<&str> = line.splitn(6, '\t').collect();
        if fields.len() < 5 {
            // Malformed line — pass through unchanged
            w.write_all(line.as_bytes())?;
            w.write_all(b"\n")?;
            return Ok(());
        }

        let chrom = fields[0];
        let pos_str = fields[1];
        // fields[2] is the ID column we'll replace
        let ref_seq = fields[3];
        let alt_field = fields[4];

        // Parse CHROM
        let chr_idx = ChrIndex::from_name(chrom);

        // Parse POS
        let pos: Option<u32> = pos_str.parse().ok();

        // Compute new ID
        let new_id = match (chr_idx, pos) {
            (Some(chr), Some(p)) if p > 0 => {
                compute_id(chr, p, ref_seq, alt_field, assembly, use_uuid)
            }
            _ => {
                // Can't encode — keep original ID
                eprintln!(
                    "Warning: cannot encode variant at {}:{} — keeping original ID",
                    chrom, pos_str
                );
                fields[2].to_string()
            }
        };

        // Reconstruct the line: CHROM \t POS \t NEW_ID \t REF \t ALT [\t rest]
        w.write_all(fields[0].as_bytes())?;
        w.write_all(b"\t")?;
        w.write_all(fields[1].as_bytes())?;
        w.write_all(b"\t")?;
        w.write_all(new_id.as_bytes())?;
        w.write_all(b"\t")?;
        w.write_all(fields[3].as_bytes())?;
        w.write_all(b"\t")?;
        // fields[4] might be just ALT, or fields[5] has the rest (QUAL onwards)
        if fields.len() == 6 {
            w.write_all(fields[4].as_bytes())?;
            w.write_all(b"\t")?;
            w.write_all(fields[5].as_bytes())?;
        } else {
            w.write_all(fields[4].as_bytes())?;
        }
        w.write_all(b"\n")?;
        Ok(())
    };

    // Process the first data line (already read during header scan)
    if let Some(ref line) = first_data_line {
        process_line(line, &mut writer)?;
        count += 1;
    }

    // Process remaining data lines
    loop {
        line_buf.clear();
        let n = reader.read_line(&mut line_buf)?;
        if n == 0 {
            break;
        }
        let trimmed = line_buf.trim_end_matches(['\n', '\r']);
        if trimmed.is_empty() {
            continue;
        }
        process_line(trimmed, &mut writer)?;
        count += 1;
    }

    // Flush and finalise
    writer.flush()?;

    // For bgzf output, we need to call finish() on the underlying bgzf writer
    // to write the EOF block.
    if let OutputWriter::Bgzf(buf_writer) = writer {
        let bgzf_writer = buf_writer
            .into_inner()
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        bgzf_writer.finish()?;
    }

    Ok(count)
}

/// Compute the ID string for a single VCF record.
///
/// Handles multi-allelic ALTs (comma-separated) — produces one UVID per
/// allele, joined with `_`. If ALT is `.`, returns `.`.
fn compute_id(
    chr: ChrIndex,
    pos: u32,
    ref_seq: &str,
    alt_field: &str,
    assembly: Assembly,
    use_uuid: bool,
) -> String {
    // ALT = "." means no alternate allele
    if alt_field == "." {
        return ".".to_string();
    }

    let alt_alleles: Vec<&str> = alt_field.split(',').collect();
    let mut ids: Vec<String> = Vec::with_capacity(alt_alleles.len());

    for alt in &alt_alleles {
        let alt_trimmed = alt.trim();
        if alt_trimmed == "." || alt_trimmed == "*" {
            // Symbolic / missing — skip encoding
            ids.push(".".to_string());
            continue;
        }

        match Uvid128::encode(
            chr,
            pos,
            ref_seq.as_bytes(),
            alt_trimmed.as_bytes(),
            assembly,
        ) {
            Some(uvid) => {
                if use_uuid {
                    ids.push(uvid.to_uuid5().to_string());
                } else {
                    ids.push(uvid.to_hex_string());
                }
            }
            None => {
                eprintln!(
                    "Warning: failed to encode UVID for {}:{} {}>{}",
                    chr, pos, ref_seq, alt_trimmed
                );
                ids.push(".".to_string());
            }
        }
    }

    ids.join("_")
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -- Assembly detection tests --

    #[test]
    fn test_detect_assembly_from_reference_grch38() {
        let headers = vec!["##reference=GRCh38".to_string()];
        assert_eq!(
            detect_assembly_from_header(&headers),
            Some(Assembly::GRCh38)
        );
    }

    #[test]
    fn test_detect_assembly_from_reference_grch37_path() {
        let headers = vec!["##reference=file:///data/ref/GRCh37.fa".to_string()];
        assert_eq!(
            detect_assembly_from_header(&headers),
            Some(Assembly::GRCh37)
        );
    }

    #[test]
    fn test_detect_assembly_from_reference_hg19() {
        let headers = vec!["##reference=hg19".to_string()];
        assert_eq!(
            detect_assembly_from_header(&headers),
            Some(Assembly::GRCh37)
        );
    }

    #[test]
    fn test_detect_assembly_from_reference_hg38_path() {
        let headers = vec!["##reference=file:///data/ref/hg38.fa".to_string()];
        assert_eq!(
            detect_assembly_from_header(&headers),
            Some(Assembly::GRCh38)
        );
    }

    #[test]
    fn test_detect_assembly_from_assembly_line() {
        let headers = vec!["##assembly=hg19".to_string()];
        assert_eq!(
            detect_assembly_from_header(&headers),
            Some(Assembly::GRCh37)
        );
    }

    #[test]
    fn test_detect_assembly_from_contig_line() {
        let headers = vec!["##contig=<ID=chr1,length=248956422,assembly=GRCh38>".to_string()];
        assert_eq!(
            detect_assembly_from_header(&headers),
            Some(Assembly::GRCh38)
        );
    }

    #[test]
    fn test_detect_assembly_none_when_missing() {
        let headers = vec![
            "##fileformat=VCFv4.3".to_string(),
            "##contig=<ID=chr1,length=248956422>".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        assert_eq!(detect_assembly_from_header(&headers), None);
    }

    #[test]
    fn test_detect_assembly_prefers_first_match() {
        let headers = vec![
            "##reference=GRCh38".to_string(),
            "##contig=<ID=chr1,length=249250621,assembly=GRCh37>".to_string(),
        ];
        // Should return GRCh38 from the first matching line
        assert_eq!(
            detect_assembly_from_header(&headers),
            Some(Assembly::GRCh38)
        );
    }

    // -- compute_id tests --

    #[test]
    fn test_compute_id_single_alt() {
        let id = compute_id(
            ChrIndex(0), // chr1
            100,
            "A",
            "G",
            Assembly::GRCh38,
            false,
        );
        // Should produce a valid UVID hex string (32 hex chars + 3 dashes = 35 chars)
        assert_eq!(id.len(), 35);
        assert!(!id.contains('_'));
    }

    #[test]
    fn test_compute_id_multiallelic() {
        let id = compute_id(
            ChrIndex(1), // chr2
            1000,
            "G",
            "A,C",
            Assembly::GRCh38,
            false,
        );
        // Should contain exactly one underscore (two alleles)
        let parts: Vec<&str> = id.split('_').collect();
        assert_eq!(parts.len(), 2);
        assert_eq!(parts[0].len(), 35);
        assert_eq!(parts[1].len(), 35);
    }

    #[test]
    fn test_compute_id_dot_alt() {
        let id = compute_id(ChrIndex(0), 100, "A", ".", Assembly::GRCh38, false);
        assert_eq!(id, ".");
    }

    #[test]
    fn test_compute_id_uuid_mode() {
        let id = compute_id(ChrIndex(0), 100, "A", "G", Assembly::GRCh38, true);
        // UUID format: 8-4-4-4-12 = 36 chars
        assert_eq!(id.len(), 36);
        // Verify it looks like a UUID
        let parts: Vec<&str> = id.split('-').collect();
        assert_eq!(parts.len(), 5);
    }

    // -- Full passthrough tests using in-memory reader --

    fn make_test_vcf(header_extra: &str) -> String {
        let mut vcf = String::new();
        vcf.push_str("##fileformat=VCFv4.3\n");
        if !header_extra.is_empty() {
            vcf.push_str(header_extra);
            if !header_extra.ends_with('\n') {
                vcf.push('\n');
            }
        }
        vcf.push_str("##contig=<ID=chr1,length=248956422>\n");
        vcf.push_str("##contig=<ID=chr2,length=242193529>\n");
        vcf.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        vcf.push_str("chr1\t100\t.\tA\tG\t30\tPASS\tDP=50\n");
        vcf.push_str("chr1\t200\trs123\tC\tT\t45\tPASS\tDP=60\n");
        vcf.push_str("chr2\t1000\t.\tG\tA,C\t50\tPASS\tDP=80\n");
        vcf.push_str("chr2\t2000\t.\tT\t.\t10\tPASS\tDP=10\n");
        vcf
    }

    #[test]
    fn test_passthrough_basic_with_assembly_override() {
        let vcf_data = make_test_vcf("");

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let out_path = tmp.path().to_path_buf();

        // Write input to a temp file
        let input_tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(input_tmp.path(), &vcf_data).unwrap();

        let count = vcf_passthrough(
            input_tmp.path(),
            Some(&out_path),
            false,
            Some(Assembly::GRCh38),
        )
        .unwrap();

        assert_eq!(count, 4);

        let output = std::fs::read_to_string(&out_path).unwrap();
        let lines: Vec<&str> = output.lines().collect();

        // Header lines should pass through unchanged
        assert!(lines[0].starts_with("##fileformat="));
        assert!(lines[1].starts_with("##contig="));
        assert!(lines[2].starts_with("##contig="));
        assert!(lines[3].starts_with("#CHROM"));

        // Data line 1: chr1 100 . A G — should have UVID
        let fields4: Vec<&str> = lines[4].split('\t').collect();
        assert_eq!(fields4[0], "chr1");
        assert_eq!(fields4[1], "100");
        assert_ne!(fields4[2], "."); // ID replaced
        assert_eq!(fields4[2].len(), 35); // UVID hex format
        assert_eq!(fields4[3], "A");
        assert_eq!(fields4[4], "G");

        // Data line 2: chr1 200 rs123 C T — original ID replaced with UVID
        let fields5: Vec<&str> = lines[5].split('\t').collect();
        assert_ne!(fields5[2], "rs123"); // replaced
        assert_eq!(fields5[2].len(), 35);

        // Data line 3: chr2 1000 . G A,C — multi-allelic, _ separated
        let fields6: Vec<&str> = lines[6].split('\t').collect();
        assert!(fields6[2].contains('_'));
        let uvid_parts: Vec<&str> = fields6[2].split('_').collect();
        assert_eq!(uvid_parts.len(), 2);

        // Data line 4: chr2 2000 . T . — ALT is ".", ID stays "."
        let fields7: Vec<&str> = lines[7].split('\t').collect();
        assert_eq!(fields7[2], ".");
    }

    #[test]
    fn test_passthrough_assembly_detection_error() {
        // VCF with no assembly info and no override → should fail
        let vcf_data = make_test_vcf("");
        let input_tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(input_tmp.path(), &vcf_data).unwrap();

        let result = vcf_passthrough(input_tmp.path(), None, false, None);
        assert!(result.is_err());
        match result.unwrap_err() {
            VcfPassthroughError::AssemblyNotDetected => {} // expected
            other => panic!("Expected AssemblyNotDetected, got: {}", other),
        }
    }

    #[test]
    fn test_passthrough_assembly_detected_from_header() {
        let vcf_data = make_test_vcf("##reference=GRCh38");
        let input_tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(input_tmp.path(), &vcf_data).unwrap();

        let out_tmp = tempfile::NamedTempFile::new().unwrap();
        let count = vcf_passthrough(
            input_tmp.path(),
            Some(out_tmp.path()),
            false,
            None, // no override — should detect from header
        )
        .unwrap();

        assert_eq!(count, 4);
    }

    #[test]
    fn test_passthrough_uuid_mode() {
        let vcf_data = make_test_vcf("##reference=GRCh38");
        let input_tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(input_tmp.path(), &vcf_data).unwrap();

        let out_tmp = tempfile::NamedTempFile::new().unwrap();
        vcf_passthrough(
            input_tmp.path(),
            Some(out_tmp.path()),
            true, // UUID mode
            None,
        )
        .unwrap();

        let output = std::fs::read_to_string(out_tmp.path()).unwrap();
        let data_lines: Vec<&str> = output.lines().filter(|l| !l.starts_with('#')).collect();

        // First record: should have UUID format (36 chars)
        let fields: Vec<&str> = data_lines[0].split('\t').collect();
        assert_eq!(fields[2].len(), 36);
        assert_eq!(fields[2].split('-').count(), 5);

        // Multi-allelic: two UUIDs joined by _
        let fields3: Vec<&str> = data_lines[2].split('\t').collect();
        let parts: Vec<&str> = fields3[2].split('_').collect();
        assert_eq!(parts.len(), 2);
        assert_eq!(parts[0].len(), 36);
        assert_eq!(parts[1].len(), 36);
    }

    #[test]
    fn test_passthrough_bgzf_output() {
        let vcf_data = make_test_vcf("##reference=GRCh38");
        let input_tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(input_tmp.path(), &vcf_data).unwrap();

        let out_dir = tempfile::tempdir().unwrap();
        let out_path = out_dir.path().join("output.vcf.gz");

        vcf_passthrough(input_tmp.path(), Some(&out_path), false, None).unwrap();

        // Verify the output file starts with gzip magic bytes
        let bytes = std::fs::read(&out_path).unwrap();
        assert!(bytes.len() > 2);
        assert_eq!(bytes[0], 0x1f);
        assert_eq!(bytes[1], 0x8b);

        // Verify we can read it back with bgzf reader
        let file = std::fs::File::open(&out_path).unwrap();
        let decoder = bgzf::io::Reader::new(file);
        let reader = BufReader::new(decoder);
        let lines: Vec<String> = reader.lines().collect::<Result<_, _>>().unwrap();
        let data_lines: Vec<&String> = lines.iter().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data_lines.len(), 4);
    }

    #[test]
    fn test_passthrough_preserves_all_columns() {
        let vcf_data = make_test_vcf("##reference=GRCh38");
        let input_tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(input_tmp.path(), &vcf_data).unwrap();

        let out_tmp = tempfile::NamedTempFile::new().unwrap();
        vcf_passthrough(input_tmp.path(), Some(out_tmp.path()), false, None).unwrap();

        let output = std::fs::read_to_string(out_tmp.path()).unwrap();
        let data_lines: Vec<&str> = output.lines().filter(|l| !l.starts_with('#')).collect();

        // First record: verify all non-ID columns are preserved
        let fields: Vec<&str> = data_lines[0].split('\t').collect();
        assert_eq!(fields[0], "chr1");
        assert_eq!(fields[1], "100");
        // fields[2] is the replaced ID
        assert_eq!(fields[3], "A");
        assert_eq!(fields[4], "G");
        assert_eq!(fields[5], "30");
        assert_eq!(fields[6], "PASS");
        assert_eq!(fields[7], "DP=50");
    }

    #[test]
    fn test_passthrough_deterministic() {
        let vcf_data = make_test_vcf("##reference=GRCh38");

        let input_tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(input_tmp.path(), &vcf_data).unwrap();

        let out1 = tempfile::NamedTempFile::new().unwrap();
        let out2 = tempfile::NamedTempFile::new().unwrap();

        vcf_passthrough(input_tmp.path(), Some(out1.path()), false, None).unwrap();
        vcf_passthrough(input_tmp.path(), Some(out2.path()), false, None).unwrap();

        let content1 = std::fs::read_to_string(out1.path()).unwrap();
        let content2 = std::fs::read_to_string(out2.path()).unwrap();
        assert_eq!(content1, content2);
    }
}
