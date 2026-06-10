//! Sample manifest: auto-detect paired FASTQs and batch process.
//!
//! Scans a directory for paired-end FASTQ files, groups them by sample name,
//! and generates a TSV manifest for reproducible batch processing.
//!
//! Supports common Illumina naming conventions:
//!   - SampleName_S1_L001_R1_001.fastq.gz / _R2_001.fastq.gz
//!   - SampleName_R1.fastq.gz / _R2.fastq.gz
//!   - SampleName_1.fastq.gz / _2.fastq.gz

use anyhow::{bail, Context, Result};
use log::{info, warn};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

/// A single sample entry in the manifest.
#[derive(Debug, Clone)]
pub struct SampleEntry {
    pub sample_name: String,
    pub fastq1: PathBuf,
    pub fastq2: PathBuf,
}

/// Scan a directory for paired FASTQ files and return sample entries.
///
/// Pairing logic:
///   1. Find all files matching *.fastq.gz, *.fq.gz, *.fastq, *.fq
///   2. Identify R1/R2 pairs by common naming patterns
///   3. Group by derived sample name
pub fn scan_directory(dir: &Path, recursive: bool) -> Result<Vec<SampleEntry>> {
    if !dir.is_dir() {
        bail!("{} is not a directory", dir.display());
    }

    // Collect all FASTQ files
    let mut fastq_files: Vec<PathBuf> = Vec::new();
    collect_fastqs(dir, &mut fastq_files, recursive)?;
    fastq_files.sort();

    if fastq_files.is_empty() {
        bail!("No FASTQ files found in {}", dir.display());
    }

    info!("Found {} FASTQ files in {}", fastq_files.len(), dir.display());

    // Group into R1/R2 pairs
    let mut r1_files: BTreeMap<String, PathBuf> = BTreeMap::new();
    let mut r2_files: BTreeMap<String, PathBuf> = BTreeMap::new();

    for path in &fastq_files {
        let fname = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("");

        if let Some((sample, read)) = parse_fastq_name(fname) {
            match read {
                ReadNumber::R1 => {
                    r1_files.insert(sample, path.clone());
                }
                ReadNumber::R2 => {
                    r2_files.insert(sample, path.clone());
                }
            }
        } else {
            warn!("Cannot determine read number for: {}", fname);
        }
    }

    // Match R1 to R2
    let mut samples = Vec::new();
    let mut unmatched = Vec::new();

    for (sample_key, r1_path) in &r1_files {
        if let Some(r2_path) = r2_files.get(sample_key) {
            // Derive a clean sample name from the key
            let sample_name = clean_sample_name(sample_key);
            samples.push(SampleEntry {
                sample_name,
                fastq1: r1_path.clone(),
                fastq2: r2_path.clone(),
            });
        } else {
            unmatched.push(r1_path.clone());
        }
    }

    // Check for R2 without R1
    for (sample_key, r2_path) in &r2_files {
        if !r1_files.contains_key(sample_key) {
            unmatched.push(r2_path.clone());
        }
    }

    if !unmatched.is_empty() {
        warn!(
            "{} unpaired FASTQ file(s) (no matching R1/R2):",
            unmatched.len()
        );
        for u in &unmatched {
            warn!("  {}", u.display());
        }
    }

    if samples.is_empty() {
        bail!(
            "No paired FASTQ samples found. Found {} R1 and {} R2 files but no matching pairs.",
            r1_files.len(),
            r2_files.len()
        );
    }

    info!("Identified {} paired sample(s)", samples.len());
    Ok(samples)
}

/// Write a manifest TSV file.
pub fn write_manifest(samples: &[SampleEntry], output: &Path) -> Result<()> {
    let mut w = BufWriter::new(File::create(output)?);
    writeln!(w, "sample_name\tfastq1\tfastq2")?;
    for s in samples {
        writeln!(
            w,
            "{}\t{}\t{}",
            s.sample_name,
            s.fastq1.display(),
            s.fastq2.display()
        )?;
    }
    w.flush()?;
    info!("Manifest written to: {}", output.display());
    Ok(())
}

/// Read a manifest TSV file.
pub fn read_manifest(path: &Path) -> Result<Vec<SampleEntry>> {
    let reader = BufReader::new(
        File::open(path).with_context(|| format!("Cannot open manifest: {}", path.display()))?,
    );

    let mut samples = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim();

        // Skip header and comments
        if i == 0 && line.starts_with("sample_name") {
            continue;
        }
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            // Also try comma-separated
            let fields: Vec<&str> = line.split(',').collect();
            if fields.len() < 3 {
                warn!("Skipping malformed line {}: {}", i + 1, line);
                continue;
            }
            samples.push(SampleEntry {
                sample_name: fields[0].trim().to_string(),
                fastq1: PathBuf::from(fields[1].trim()),
                fastq2: PathBuf::from(fields[2].trim()),
            });
            continue;
        }

        samples.push(SampleEntry {
            sample_name: fields[0].trim().to_string(),
            fastq1: PathBuf::from(fields[1].trim()),
            fastq2: PathBuf::from(fields[2].trim()),
        });
    }

    if samples.is_empty() {
        bail!("No samples found in manifest: {}", path.display());
    }

    info!("Loaded {} sample(s) from manifest", samples.len());
    Ok(samples)
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

#[derive(Debug, PartialEq)]
enum ReadNumber {
    R1,
    R2,
}

/// Parse an Illumina-style FASTQ filename into (sample_key, read_number).
///
/// Supported patterns:
///   "SampleName_S1_L001_R1_001.fastq.gz" -> ("SampleName_S1_L001", R1)
///   "SampleName_R1_001.fastq.gz"         -> ("SampleName", R1)
///   "SampleName_R1.fastq.gz"             -> ("SampleName", R1)
///   "SampleName_1.fastq.gz"              -> ("SampleName", R1)
fn parse_fastq_name(filename: &str) -> Option<(String, ReadNumber)> {
    // Strip extensions: .fastq.gz, .fq.gz, .fastq, .fq
    let base = filename
        .strip_suffix(".fastq.gz")
        .or_else(|| filename.strip_suffix(".fq.gz"))
        .or_else(|| filename.strip_suffix(".fastq"))
        .or_else(|| filename.strip_suffix(".fq"))?;

    // Try patterns in order of specificity

    // Pattern 1: _R1_001 / _R2_001 (Illumina standard)
    if let Some(prefix) = base.strip_suffix("_R1_001") {
        return Some((prefix.to_string(), ReadNumber::R1));
    }
    if let Some(prefix) = base.strip_suffix("_R2_001") {
        return Some((prefix.to_string(), ReadNumber::R2));
    }

    // Pattern 2: _R1 / _R2 (common shorthand)
    if let Some(prefix) = base.strip_suffix("_R1") {
        return Some((prefix.to_string(), ReadNumber::R1));
    }
    if let Some(prefix) = base.strip_suffix("_R2") {
        return Some((prefix.to_string(), ReadNumber::R2));
    }

    // Pattern 3: _1 / _2 (minimal)
    if let Some(prefix) = base.strip_suffix("_1") {
        return Some((prefix.to_string(), ReadNumber::R1));
    }
    if let Some(prefix) = base.strip_suffix("_2") {
        return Some((prefix.to_string(), ReadNumber::R2));
    }

    // Pattern 4: .R1 / .R2 (dot-separated)
    if let Some(prefix) = base.strip_suffix(".R1") {
        return Some((prefix.to_string(), ReadNumber::R1));
    }
    if let Some(prefix) = base.strip_suffix(".R2") {
        return Some((prefix.to_string(), ReadNumber::R2));
    }

    None
}

/// Derive a clean sample name from the pairing key.
///
/// "25NGS1071-Duplex_S1_L001" -> "25NGS1071-Duplex"
/// Strips _S\d+, _L\d+ Illumina suffixes for cleaner naming.
fn clean_sample_name(key: &str) -> String {
    let mut name = key.to_string();

    // Strip trailing _L00N lane identifiers
    if let Some(idx) = name.rfind("_L00") {
        if name[idx..].len() >= 5
            && name[idx + 4..].chars().next().map_or(false, |c| c.is_ascii_digit())
        {
            name.truncate(idx);
        }
    }

    // Strip trailing _SN sample index
    if let Some(idx) = name.rfind("_S") {
        let suffix = &name[idx + 2..];
        if !suffix.is_empty() && suffix.chars().all(|c| c.is_ascii_digit()) {
            name.truncate(idx);
        }
    }

    // If stripping left us empty, use the original key
    if name.is_empty() {
        name = key.to_string();
    }

    name
}

/// Recursively collect FASTQ files from a directory.
fn collect_fastqs(dir: &Path, out: &mut Vec<PathBuf>, recursive: bool) -> Result<()> {
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();

        if path.is_dir() && recursive {
            collect_fastqs(&path, out, recursive)?;
        } else if path.is_file() {
            let name = path
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or("");
            if name.ends_with(".fastq.gz")
                || name.ends_with(".fq.gz")
                || name.ends_with(".fastq")
                || name.ends_with(".fq")
            {
                out.push(path);
            }
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_illumina_standard() {
        let (key, read) =
            parse_fastq_name("25NGS1071-Duplex_S1_L001_R1_001.fastq.gz").unwrap();
        assert_eq!(key, "25NGS1071-Duplex_S1_L001");
        assert_eq!(read, ReadNumber::R1);

        let (key, read) =
            parse_fastq_name("25NGS1071-Duplex_S1_L001_R2_001.fastq.gz").unwrap();
        assert_eq!(key, "25NGS1071-Duplex_S1_L001");
        assert_eq!(read, ReadNumber::R2);
    }

    #[test]
    fn test_parse_short_names() {
        let (key, read) = parse_fastq_name("sampleA_R1.fastq.gz").unwrap();
        assert_eq!(key, "sampleA");
        assert_eq!(read, ReadNumber::R1);

        let (key, read) = parse_fastq_name("sampleA_1.fq.gz").unwrap();
        assert_eq!(key, "sampleA");
        assert_eq!(read, ReadNumber::R1);
    }

    #[test]
    fn test_clean_sample_name() {
        assert_eq!(
            clean_sample_name("25NGS1071-Duplex_S1_L001"),
            "25NGS1071-Duplex"
        );
        assert_eq!(clean_sample_name("SampleA_S12"), "SampleA");
        assert_eq!(clean_sample_name("MySample"), "MySample");
    }

    #[test]
    fn test_unparseable() {
        assert!(parse_fastq_name("readme.txt").is_none());
        assert!(parse_fastq_name("sample.bam").is_none());
    }
}
