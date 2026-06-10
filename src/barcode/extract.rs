//! Extract UMI barcodes from paired-end FASTQ reads and align with configurable backend.
//!
//! Supported aligner backends:
//!   - **Parabricks** (GPU): `pbrun fq2bam` -- GPU-accelerated BWA-MEM + sort + dedup
//!   - **BWA-MEM2** (CPU):   SIMD-optimized BWA-MEM, 2-3x faster than BWA
//!   - **minimap2** (CPU):   3x faster than BWA for short reads, standard for long reads
//!   - **bwa** (CPU):        Legacy fallback
//!
//! Workflow:
//!   1. Parse barcode pattern to determine UMI vs spacer positions
//!   2. Read paired FASTQ files
//!   3. Extract UMI bases, validate spacer bases
//!   4. Concatenate R1+R2 UMIs, append to read header
//!   5. Trim UMI+spacer from read sequence and quality
//!   6. Write trimmed FASTQs
//!   7. Invoke chosen aligner -> sorted, indexed BAM

use anyhow::{Context, Result, bail};
use log::{info, warn};
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::str::FromStr;

// ---------------------------------------------------------------------------
// Aligner backend configuration
// ---------------------------------------------------------------------------

/// Which aligner to use for the FASTQ -> sorted BAM step.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignerBackend {
    /// NVIDIA Clara Parabricks `pbrun fq2bam` -- full GPU acceleration.
    /// Replaces BWA-MEM + samtools sort + MarkDuplicates in a single GPU pass.
    /// Free for research use; requires NVIDIA GPU with >= 16 GB VRAM.
    Parabricks,

    /// BWA-MEM2 -- SIMD-optimized CPU rewrite of BWA-MEM.
    /// Drop-in compatible output, 2-3x faster than original BWA.
    /// Open source (MIT), no GPU required.
    BwaMem2,

    /// minimap2 -- fast general-purpose aligner.
    /// 3x faster than BWA-MEM for Illumina short reads.
    /// Standard for long reads (PacBio, ONT). Open source (MIT).
    Minimap2,

    /// Original BWA-MEM -- legacy fallback.
    /// Slowest option; a warning is printed recommending alternatives.
    Bwa,
}

impl fmt::Display for AlignerBackend {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parabricks => write!(f, "parabricks"),
            Self::BwaMem2 => write!(f, "bwa-mem2"),
            Self::Minimap2 => write!(f, "minimap2"),
            Self::Bwa => write!(f, "bwa"),
        }
    }
}

impl FromStr for AlignerBackend {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        match s.to_lowercase().as_str() {
            "parabricks" | "pbrun" | "pb" => Ok(Self::Parabricks),
            "bwa-mem2" | "bwamem2" => Ok(Self::BwaMem2),
            "minimap2" | "mm2" => Ok(Self::Minimap2),
            "bwa" | "bwa-mem" => Ok(Self::Bwa),
            _ => bail!(
                "Unknown aligner '{}'. Choose from: parabricks, bwa-mem2, minimap2, bwa",
                s
            ),
        }
    }
}

impl AlignerBackend {
    /// Default executable name for this backend.
    pub fn default_executable(&self) -> &'static str {
        match self {
            Self::Parabricks => "pbrun",
            Self::BwaMem2 => "bwa-mem2",
            Self::Minimap2 => "minimap2",
            Self::Bwa => "bwa",
        }
    }

    /// Human-readable description for log messages.
    pub fn description(&self) -> &'static str {
        match self {
            Self::Parabricks => "NVIDIA Clara Parabricks (GPU-accelerated BWA-MEM)",
            Self::BwaMem2 => "BWA-MEM2 (SIMD-optimized CPU)",
            Self::Minimap2 => "minimap2 (fast CPU aligner)",
            Self::Bwa => "BWA-MEM (legacy CPU)",
        }
    }
}

/// Full configuration for the alignment step.
#[derive(Debug, Clone)]
pub struct AlignerConfig {
    pub backend: AlignerBackend,
    /// Override path to aligner executable (None = use default from PATH).
    pub executable_path: Option<String>,
    /// Reference genome FASTA (must be indexed for chosen aligner).
    pub reference: PathBuf,
    /// Extra arguments passed through to the aligner.
    pub extra_args: Option<String>,
    /// Use Docker to invoke Parabricks (instead of native pbrun).
    pub parabricks_docker: bool,
    /// Parabricks Docker image tag.
    pub parabricks_image: String,
    /// Path to samtools binary (for sort/index with non-Parabricks backends).
    pub samtools_path: String,
    /// Number of CPU threads.
    pub threads: usize,
    /// CUDA device index (Parabricks only).
    pub gpu_device: usize,
}

impl AlignerConfig {
    /// Resolve the actual executable path.
    fn executable(&self) -> &str {
        self.executable_path
            .as_deref()
            .unwrap_or(self.backend.default_executable())
    }

    /// Check if the executable is reachable.
    ///
    /// Uses `which` to locate the binary in PATH. We cannot rely on
    /// running `exe --version` because some bioinformatics tools
    /// (notably bwa-mem2) return non-zero exit codes for --version.
    fn check_executable(&self) -> Result<()> {
        let exe = self.executable();

        // For Docker-mode Parabricks, check Docker instead of pbrun.
        if self.backend == AlignerBackend::Parabricks && self.parabricks_docker {
            let status = Command::new("which")
                .arg("docker")
                .stdout(std::process::Stdio::null())
                .stderr(std::process::Stdio::null())
                .status();
            match status {
                Ok(s) if s.success() => return Ok(()),
                _ => bail!(
                    "Docker is required for Parabricks Docker mode but was not found. \
                     Install Docker or use native pbrun with --aligner-path."
                ),
            }
        }

        // First try `which exe` to check PATH
        let which_result = Command::new("which")
            .arg(exe)
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status();

        // Also accept an absolute/relative path that exists on disk
        let found = match which_result {
            Ok(s) if s.success() => true,
            _ => std::path::Path::new(exe).exists(),
        };

        if found {
            Ok(())
        } else {
            let hint = match self.backend {
                AlignerBackend::Parabricks => {
                    "Install from: https://docs.nvidia.com/clara/parabricks/latest/\n\
                     Or use Docker: --parabricks-docker"
                }
                AlignerBackend::BwaMem2 => {
                    "Install: conda install -c bioconda bwa-mem2\n\
                     Or: https://github.com/bwa-mem2/bwa-mem2"
                }
                AlignerBackend::Minimap2 => {
                    "Install: conda install -c bioconda minimap2\n\
                     Or: https://github.com/lh3/minimap2"
                }
                AlignerBackend::Bwa => {
                    "Install: conda install -c bioconda bwa\n\
                     Consider using --aligner bwa-mem2 for 2-3x faster alignment."
                }
            };
            bail!(
                "Aligner '{}' not found at '{}'. {}\n{}",
                self.backend,
                exe,
                "Ensure it is installed and in your PATH, or set --aligner-path.",
                hint
            )
        }
    }
}

// ---------------------------------------------------------------------------
// Barcode pattern parsing
// ---------------------------------------------------------------------------

/// Parsed barcode pattern: each position is UMI, validated spacer, or skip.
///
/// Supports IDT fgbio-style read structures:
///   N = random UMI base (extracted into barcode)
///   S = skip base (trimmed from read, not validated, not part of UMI)
///   A/C/G/T = fixed spacer base (trimmed, validated against expected sequence)
///
/// Example: "NNNSS" for IDT xGen Duplex Seq (3bp UMI + 2bp skip).
#[derive(Debug, Clone)]
pub struct BarcodePattern {
    pub positions: Vec<bool>,
    pub expected_bases: Vec<Option<u8>>,
    pub total_len: usize,
    pub umi_len: usize,
}

impl BarcodePattern {
    /// Parse a barcode pattern string.
    /// N = random UMI base, S = skip (trim but don't validate), A/C/G/T = fixed spacer.
    /// Example: "NNNSS" for IDT xGen Duplex Seq (3bp UMI + 2bp skip per read).
    /// Example: "NNNGT" means 3bp UMI followed by validated "GT" spacer.
    pub fn from_pattern(pattern: &str) -> Result<Self> {
        let mut positions = Vec::with_capacity(pattern.len());
        let mut expected_bases = Vec::with_capacity(pattern.len());
        let mut umi_len = 0;

        for ch in pattern.chars() {
            match ch {
                'N' | 'n' => {
                    positions.push(true);
                    expected_bases.push(None);
                    umi_len += 1;
                }
                'S' | 's' => {
                    // Skip base: trimmed from read but not validated or included in UMI.
                    // Matches IDT fgbio read structure notation (e.g., 3M2S146T).
                    positions.push(false);
                    expected_bases.push(None);
                }
                'A' | 'C' | 'G' | 'T' => {
                    positions.push(false);
                    expected_bases.push(Some(ch as u8));
                }
                _ => bail!(
                    "Invalid barcode pattern character '{}'. Use N for UMI, S for skip, A/C/G/T for fixed spacer.",
                    ch
                ),
            }
        }

        Ok(Self {
            total_len: pattern.len(),
            positions,
            expected_bases,
            umi_len,
        })
    }

    /// Extract UMI bases from a sequence, validate spacer positions.
    pub fn extract_umi(&self, seq: &[u8]) -> (String, bool) {
        if seq.len() < self.total_len {
            return (String::new(), false);
        }

        let mut umi = Vec::with_capacity(self.umi_len);
        let mut spacer_valid = true;

        for (i, (&is_umi, &expected)) in
            self.positions.iter().zip(self.expected_bases.iter()).enumerate()
        {
            if is_umi {
                umi.push(seq[i]);
            } else if let Some(exp) = expected {
                if seq[i] != exp {
                    spacer_valid = false;
                }
            }
        }

        (String::from_utf8_lossy(&umi).to_string(), spacer_valid)
    }
}

/// Statistics from barcode extraction.
#[derive(Debug, Default)]
pub struct ExtractionStats {
    pub total_pairs: u64,
    pub passing_pairs: u64,
    pub bad_spacer: u64,
    pub bad_barcode: u64,
}

// ---------------------------------------------------------------------------
// Main pipeline entry point
// ---------------------------------------------------------------------------

/// Run the full fastq-to-bam pipeline:
///   1. Extract barcodes from paired FASTQs
///   2. Align with chosen backend
///   3. Sort and index BAM
pub fn run_fastq2bam(
    fastq1: &Path,
    fastq2: &Path,
    output_dir: &Path,
    bpattern: Option<&str>,
    blist: Option<&Path>,
    config: &AlignerConfig,
) -> Result<()> {
    std::fs::create_dir_all(output_dir)?;

    // Validate aligner availability
    config.check_executable()?;

    // Print performance advisory for legacy BWA users
    if config.backend == AlignerBackend::Bwa {
        warn!(
            "Using legacy BWA-MEM. Consider faster alternatives:\n  \
             --aligner bwa-mem2    : 2-3x faster (CPU, drop-in replacement)\n  \
             --aligner minimap2    : 3x faster for Illumina short reads (CPU)\n  \
             --aligner parabricks  : 10-60x faster (GPU, NVIDIA CUDA required)"
        );
    }

    info!("Aligner: {}", config.backend.description());

    // Parse barcode specification
    let pattern = match (bpattern, blist) {
        (Some(pat), _) => BarcodePattern::from_pattern(pat)?,
        (None, Some(_blist_path)) => {
            info!("Barcode list mode: loading valid barcodes");
            // TODO: implement barcode list validation
            BarcodePattern::from_pattern("NNN")?
        }
        (None, None) => bail!("Must specify either --bpattern or --blist"),
    };

    info!(
        "Barcode pattern: {} total bases, {} UMI bases, {} spacer bases",
        pattern.total_len,
        pattern.umi_len,
        pattern.total_len - pattern.umi_len
    );

    // Derive output paths
    let stem = fastq1
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("sample")
        .trim_end_matches(".fastq")
        .trim_end_matches(".fq");

    let tagged_fq1 = output_dir.join(format!("{}_tagged_R1.fastq", stem));
    let tagged_fq2 = output_dir.join(format!("{}_tagged_R2.fastq", stem));
    let out_bam = output_dir.join(format!("{}.bam", stem));
    let stats_file = output_dir.join(format!("{}_barcode_stats.txt", stem));

    // Stage 1: Extract barcodes
    info!("Extracting barcodes from FASTQ pair...");
    let stats = extract_barcodes_from_fastq(fastq1, fastq2, &tagged_fq1, &tagged_fq2, &pattern)?;

    info!(
        "Extraction complete: {}/{} pairs passed ({:.1}%), {} bad spacer, {} bad barcode",
        stats.passing_pairs,
        stats.total_pairs,
        100.0 * stats.passing_pairs as f64 / stats.total_pairs.max(1) as f64,
        stats.bad_spacer,
        stats.bad_barcode,
    );

    let mut sf = File::create(&stats_file)?;
    writeln!(sf, "Total pairs: {}", stats.total_pairs)?;
    writeln!(sf, "Passing pairs: {}", stats.passing_pairs)?;
    writeln!(sf, "Bad spacer: {}", stats.bad_spacer)?;
    writeln!(sf, "Bad barcode: {}", stats.bad_barcode)?;

    // Stage 2: Align
    info!("Aligning with {}...", config.backend);
    match config.backend {
        AlignerBackend::Parabricks => {
            run_parabricks_alignment(&tagged_fq1, &tagged_fq2, &out_bam, config)?;
        }
        AlignerBackend::BwaMem2 | AlignerBackend::Bwa => {
            run_bwamem_style_alignment(&tagged_fq1, &tagged_fq2, &out_bam, config)?;
        }
        AlignerBackend::Minimap2 => {
            run_minimap2_alignment(&tagged_fq1, &tagged_fq2, &out_bam, config)?;
        }
    }

    // Clean up intermediate tagged FASTQs (typically 30-60 GB each)
    // These are no longer needed once alignment produces the BAM.
    if out_bam.exists() {
        info!("Alignment successful, cleaning up intermediate tagged FASTQs...");
        for fq in &[&tagged_fq1, &tagged_fq2] {
            match std::fs::remove_file(fq) {
                Ok(_) => info!("  Removed: {}", fq.display()),
                Err(e) => log::warn!("  Could not remove {}: {}", fq.display(), e),
            }
        }
    } else {
        log::warn!("Alignment output not found, keeping intermediate FASTQs for debugging");
    }

    info!("Output BAM: {}", out_bam.display());
    Ok(())
}

// ---------------------------------------------------------------------------
// Aligner backend implementations
// ---------------------------------------------------------------------------

/// Run NVIDIA Clara Parabricks `pbrun fq2bam`.
///
/// Parabricks performs GPU-accelerated:
///   - BWA-MEM alignment
///   - Coordinate sorting
///   - Duplicate marking (disabled here; UMI families handle dedup)
///   - BAM indexing
///
/// All in a single pass on the GPU -- no samtools needed.
fn run_parabricks_alignment(
    fq1: &Path,
    fq2: &Path,
    out_bam: &Path,
    config: &AlignerConfig,
) -> Result<()> {
    let out_dir = out_bam.parent().unwrap_or(Path::new("."));

    if config.parabricks_docker {
        info!("Using Parabricks Docker image: {}", config.parabricks_image);

        let workdir = out_dir.canonicalize()?;
        let refdir = config
            .reference
            .parent()
            .unwrap_or(Path::new("."))
            .canonicalize()?;
        let ref_name = config.reference.file_name().unwrap().to_str().unwrap();

        let mut cmd = Command::new("docker");
        cmd.args([
            "run",
            "--rm",
            "--gpus",
            "all",
            "--volume",
            &format!("{}:/workdir", workdir.display()),
            "--volume",
            &format!("{}:/refdir", refdir.display()),
            "--workdir",
            "/workdir",
            &config.parabricks_image,
            "pbrun",
            "fq2bam",
            "--ref",
            &format!("/refdir/{}", ref_name),
            "--in-fq",
            &format!("/workdir/{}", fq1.file_name().unwrap().to_str().unwrap()),
            &format!("/workdir/{}", fq2.file_name().unwrap().to_str().unwrap()),
            "--out-bam",
            &format!("/workdir/{}", out_bam.file_name().unwrap().to_str().unwrap()),
            "--no-markdups",
        ]);

        if let Some(ref extra) = config.extra_args {
            for arg in extra.split_whitespace() {
                cmd.arg(arg);
            }
        }

        let status = cmd
            .status()
            .context("Failed to run Parabricks via Docker")?;
        if !status.success() {
            bail!(
                "Parabricks Docker fq2bam failed with exit code {:?}",
                status.code()
            );
        }
    } else {
        let exe = config.executable();
        let mut cmd = Command::new(exe);
        cmd.args(["fq2bam", "--ref"]);
        cmd.arg(config.reference.to_str().unwrap());
        cmd.arg("--in-fq");
        cmd.arg(fq1.to_str().unwrap());
        cmd.arg(fq2.to_str().unwrap());
        cmd.arg("--out-bam");
        cmd.arg(out_bam.to_str().unwrap());
        cmd.arg("--no-markdups");
        cmd.args(["--num-gpus", "1"]);
        cmd.args(["--gpu-devices", &config.gpu_device.to_string()]);

        if let Some(ref extra) = config.extra_args {
            for arg in extra.split_whitespace() {
                cmd.arg(arg);
            }
        }

        info!("Running: {:?}", cmd);
        let status = cmd.status().context("Failed to run pbrun fq2bam")?;
        if !status.success() {
            bail!(
                "pbrun fq2bam failed with exit code {:?}",
                status.code()
            );
        }
    }

    // Parabricks produces sorted + indexed BAM directly; verify index.
    let bai_path = format!("{}.bai", out_bam.display());
    if !Path::new(&bai_path).exists() {
        info!("BAM index not found, generating with samtools...");
        run_samtools_index(out_bam, &config.samtools_path)?;
    }

    Ok(())
}

/// Run BWA-MEM or BWA-MEM2 alignment.
///
/// Both use identical CLI syntax. Output is piped through `samtools sort`.
fn run_bwamem_style_alignment(
    fq1: &Path,
    fq2: &Path,
    out_bam: &Path,
    config: &AlignerConfig,
) -> Result<()> {
    let exe = config.executable();
    let threads_str = config.threads.to_string();

    let mut cmd = Command::new(exe);
    cmd.args(["mem", "-t", &threads_str]);
    cmd.arg("-M"); // Mark shorter split hits as secondary
    cmd.arg("-C"); // Append FASTQ comment to SAM (preserves UMI in header)

    // Set chunk size for reproducibility with BWA-MEM2
    if config.backend == AlignerBackend::BwaMem2 {
        cmd.args(["-K", "10000000"]);
    }

    if let Some(ref extra) = config.extra_args {
        for arg in extra.split_whitespace() {
            cmd.arg(arg);
        }
    }

    cmd.arg(config.reference.to_str().unwrap());
    cmd.arg(fq1.to_str().unwrap());
    cmd.arg(fq2.to_str().unwrap());
    cmd.stdout(std::process::Stdio::piped());

    info!("Running: {} mem ...", exe);
    let aligner = cmd
        .spawn()
        .with_context(|| format!("Failed to start {}", exe))?;

    // Pipe through samtools sort
    let sort_threads = (config.threads / 4).max(1).to_string();
    let sort_status = Command::new(&config.samtools_path)
        .args(["sort", "-@", &sort_threads, "-o"])
        .arg(out_bam.to_str().unwrap())
        .stdin(aligner.stdout.unwrap())
        .status()
        .context("Failed to run samtools sort")?;

    if !sort_status.success() {
        bail!("samtools sort failed");
    }

    run_samtools_index(out_bam, &config.samtools_path)?;
    Ok(())
}

/// Run minimap2 alignment.
///
/// minimap2 CLI differs from BWA:
///   minimap2 -a -x sr --MD -t THREADS ref.fa R1.fq R2.fq | samtools sort
///
/// The `-x sr` preset is for Illumina short reads.
fn run_minimap2_alignment(
    fq1: &Path,
    fq2: &Path,
    out_bam: &Path,
    config: &AlignerConfig,
) -> Result<()> {
    let exe = config.executable();
    let threads_str = config.threads.to_string();

    let mut cmd = Command::new(exe);
    cmd.args([
        "-a",           // SAM output
        "-x", "sr",     // Short-read preset (Illumina)
        "--MD",         // MD tag for downstream tools
        "-t", &threads_str,
        "--secondary=no",
    ]);

    if let Some(ref extra) = config.extra_args {
        for arg in extra.split_whitespace() {
            cmd.arg(arg);
        }
    }

    cmd.arg(config.reference.to_str().unwrap());
    cmd.arg(fq1.to_str().unwrap());
    cmd.arg(fq2.to_str().unwrap());
    cmd.stdout(std::process::Stdio::piped());

    info!("Running: minimap2 -a -x sr ...");
    let aligner = cmd
        .spawn()
        .with_context(|| format!("Failed to start {}", exe))?;

    let sort_threads = (config.threads / 4).max(1).to_string();
    let sort_status = Command::new(&config.samtools_path)
        .args(["sort", "-@", &sort_threads, "-o"])
        .arg(out_bam.to_str().unwrap())
        .stdin(aligner.stdout.unwrap())
        .status()
        .context("Failed to run samtools sort")?;

    if !sort_status.success() {
        bail!("samtools sort failed");
    }

    run_samtools_index(out_bam, &config.samtools_path)?;
    Ok(())
}

/// Generate BAM index with samtools.
fn run_samtools_index(bam: &Path, samtools_path: &str) -> Result<()> {
    let status = Command::new(samtools_path)
        .args(["index"])
        .arg(bam.to_str().unwrap())
        .status()
        .context("Failed to run samtools index")?;

    if !status.success() {
        bail!("samtools index failed");
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Barcode extraction (shared across all backends)
// ---------------------------------------------------------------------------

fn extract_barcodes_from_fastq(
    fq1_path: &Path,
    fq2_path: &Path,
    out1_path: &Path,
    out2_path: &Path,
    pattern: &BarcodePattern,
) -> Result<ExtractionStats> {
    let reader1 = open_fastq_reader(fq1_path)?;
    let reader2 = open_fastq_reader(fq2_path)?;
    let mut writer1 = BufWriter::new(File::create(out1_path)?);
    let mut writer2 = BufWriter::new(File::create(out2_path)?);

    let mut stats = ExtractionStats::default();
    let mut lines1 = reader1.lines();
    let mut lines2 = reader2.lines();

    loop {
        let r1 = match read_fastq_record(&mut lines1) {
            Some(r) => r,
            None => break,
        };
        let r2 = match read_fastq_record(&mut lines2) {
            Some(r) => r,
            None => {
                warn!("Unequal number of reads in paired FASTQs");
                break;
            }
        };

        stats.total_pairs += 1;

        let (umi1, spacer1_ok) = pattern.extract_umi(r1.seq.as_bytes());
        let (umi2, spacer2_ok) = pattern.extract_umi(r2.seq.as_bytes());

        if !spacer1_ok || !spacer2_ok {
            stats.bad_spacer += 1;
            continue;
        }

        if umi1.is_empty() || umi2.is_empty() {
            stats.bad_barcode += 1;
            continue;
        }

        let combined_umi = format!("{}+{}", umi1, umi2);

        let trimmed_seq1 = &r1.seq[pattern.total_len..];
        let trimmed_qual1 = &r1.qual[pattern.total_len..];
        let trimmed_seq2 = &r2.seq[pattern.total_len..];
        let trimmed_qual2 = &r2.qual[pattern.total_len..];

        let name1_trimmed = r1.name.split_whitespace().next().unwrap_or(&r1.name);
        let name2_trimmed = r2.name.split_whitespace().next().unwrap_or(&r2.name);

        writeln!(writer1, "{}|{}", name1_trimmed, combined_umi)?;
        writeln!(writer1, "{}", trimmed_seq1)?;
        writeln!(writer1, "+")?;
        writeln!(writer1, "{}", trimmed_qual1)?;

        writeln!(writer2, "{}|{}", name2_trimmed, combined_umi)?;
        writeln!(writer2, "{}", trimmed_seq2)?;
        writeln!(writer2, "+")?;
        writeln!(writer2, "{}", trimmed_qual2)?;

        stats.passing_pairs += 1;
    }

    Ok(stats)
}

struct FastqRecord {
    name: String,
    seq: String,
    qual: String,
}

fn read_fastq_record(
    lines: &mut impl Iterator<Item = std::io::Result<String>>,
) -> Option<FastqRecord> {
    let name = lines.next()?.ok()?;
    let seq = lines.next()?.ok()?;
    let _plus = lines.next()?.ok()?;
    let qual = lines.next()?.ok()?;
    Some(FastqRecord { name, seq, qual })
}

fn open_fastq_reader(path: &Path) -> Result<BufReader<Box<dyn std::io::Read>>> {
    let file = File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;

    let reader: Box<dyn std::io::Read> = if path
        .extension()
        .map_or(false, |ext| ext == "gz" || ext == "gzip")
    {
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };

    Ok(BufReader::with_capacity(1024 * 1024, reader))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_barcode_pattern_nnn() {
        let pat = BarcodePattern::from_pattern("NNN").unwrap();
        assert_eq!(pat.umi_len, 3);
        assert_eq!(pat.total_len, 3);

        let (umi, valid) = pat.extract_umi(b"ACGTTTTT");
        assert_eq!(umi, "ACG");
        assert!(valid);
    }

    #[test]
    fn test_barcode_pattern_with_spacer() {
        let pat = BarcodePattern::from_pattern("NNNGT").unwrap();
        assert_eq!(pat.umi_len, 3);
        assert_eq!(pat.total_len, 5);

        let (umi, valid) = pat.extract_umi(b"ACGGTAAAAA");
        assert_eq!(umi, "ACG");
        assert!(valid);

        let (_, valid) = pat.extract_umi(b"ACGTTAAAAA");
        assert!(!valid);
    }

    #[test]
    fn test_aligner_backend_parsing() {
        assert_eq!(
            "parabricks".parse::<AlignerBackend>().unwrap(),
            AlignerBackend::Parabricks
        );
        assert_eq!(
            "bwa-mem2".parse::<AlignerBackend>().unwrap(),
            AlignerBackend::BwaMem2
        );
        assert_eq!(
            "minimap2".parse::<AlignerBackend>().unwrap(),
            AlignerBackend::Minimap2
        );
        assert_eq!(
            "mm2".parse::<AlignerBackend>().unwrap(),
            AlignerBackend::Minimap2
        );
        assert_eq!(
            "bwa".parse::<AlignerBackend>().unwrap(),
            AlignerBackend::Bwa
        );
        assert!("invalid".parse::<AlignerBackend>().is_err());
    }

    #[test]
    fn test_barcode_pattern_with_skip() {
        let pat = BarcodePattern::from_pattern("NNNSS").unwrap();
        assert_eq!(pat.umi_len, 3);
        assert_eq!(pat.total_len, 5);

        // Skip bases can be anything -- no validation
        let (umi, valid) = pat.extract_umi(b"ACGTTAAAAA");
        assert_eq!(umi, "ACG");
        assert!(valid);

        let (umi, valid) = pat.extract_umi(b"ACGAAAAAAA");
        assert_eq!(umi, "ACG");
        assert!(valid);
    }

    #[test]
    fn test_default_executables() {
        assert_eq!(AlignerBackend::Parabricks.default_executable(), "pbrun");
        assert_eq!(AlignerBackend::BwaMem2.default_executable(), "bwa-mem2");
        assert_eq!(AlignerBackend::Minimap2.default_executable(), "minimap2");
        assert_eq!(AlignerBackend::Bwa.default_executable(), "bwa");
    }
}
