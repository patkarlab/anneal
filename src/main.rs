//! Anneal v2: GPU-accelerated duplex consensus sequence generation.
//!
//! A production-grade Rust reimplementation of the ConsensusCruncher pipeline
//! with CUDA GPU acceleration for 10x+ speedup over the original Python.
//!
//! Pipeline stages:
//!   1. Barcode extraction from FASTQ (IDT xGen Duplex Seq adapters)
//!   2. Alignment (Parabricks GPU / BWA-MEM2 / minimap2 / BWA)
//!   3. SSCS generation (Single-Strand Consensus Sequences)
//!   4. Singleton correction (rescue unpaired reads)
//!   5. DCS generation (Duplex Consensus Sequences)

mod barcode;
mod consensus;
mod grouping;
mod manifest;
mod singleton;
mod utils;

#[cfg(feature = "gpu")]
mod cuda;

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use log::{info, warn};
use std::path::{Path, PathBuf};

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

#[derive(Parser)]
#[command(
    name = "anneal",
    version,
    about = "GPU-accelerated duplex consensus sequencing",
    long_about = "Anneal: a production-grade pipeline for UMI-based error suppression \
                  in duplex sequencing. Extracts barcodes, aligns reads, and generates \
                  high-fidelity duplex consensus sequences with optional GPU acceleration."
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    /// Number of CPU threads (default: all available)
    #[arg(long, global = true)]
    threads: Option<usize>,

    /// CUDA device index
    #[arg(long, global = true, default_value = "0")]
    gpu_device: usize,

    /// Disable GPU acceleration (CPU-only mode)
    #[arg(long, global = true)]
    no_gpu: bool,

    /// Verbosity level (-v, -vv, -vvv)
    #[arg(short, long, global = true, action = clap::ArgAction::Count)]
    verbose: u8,
}

#[derive(Subcommand)]
enum Commands {
    /// Full pipeline: extract barcodes, align, and generate consensus
    Run {
        /// Read 1 FASTQ (gzipped)
        #[arg(long)]
        fastq1: PathBuf,

        /// Read 2 FASTQ (gzipped)
        #[arg(long)]
        fastq2: PathBuf,

        /// Output directory (stages write to subdirectories)
        #[arg(short, long)]
        output: PathBuf,

        /// Barcode pattern (N = UMI, S = skip)
        /// IDT xGen Duplex Seq: "NNNSS" (3bp UMI + 2bp spacer)
        #[arg(long, default_value = "NNNSS")]
        bpattern: String,

        /// File containing list of valid barcodes
        #[arg(long)]
        blist: Option<PathBuf>,

        /// Aligner backend: "parabricks", "bwa-mem2", "minimap2", "bwa"
        #[arg(long, default_value = "bwa-mem2", value_parser = ["parabricks", "bwa-mem2", "minimap2", "bwa"])]
        aligner: String,

        /// Path to aligner executable (auto-detected if not set)
        #[arg(long)]
        aligner_path: Option<String>,

        /// Reference genome (FASTA, indexed)
        #[arg(long)]
        reference: PathBuf,

        /// Extra arguments passed to the aligner
        #[arg(long)]
        aligner_args: Option<String>,

        /// For Parabricks: use Docker instead of native pbrun
        #[arg(long)]
        parabricks_docker: bool,

        /// For Parabricks: Docker image tag
        #[arg(long, default_value = "nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1")]
        parabricks_image: String,

        /// Path to samtools
        #[arg(long, default_value = "samtools")]
        samtools: String,

        /// BED file for target regions (optional; processes all if omitted)
        #[arg(long)]
        bedfile: Option<PathBuf>,

        /// Consensus cutoff fraction (0.0-1.0)
        #[arg(long, default_value = "0.7")]
        cutoff: f64,

        /// Minimum base quality for consensus
        #[arg(long, default_value = "30")]
        min_qual: u8,

        /// Minimum reads per family for SSCS
        #[arg(long, default_value = "1")]
        min_family_size: usize,

        /// Maximum family size (downsample beyond this)
        #[arg(long, default_value = "1000")]
        max_family_size: usize,

        /// Enable singleton correction
        #[arg(long)]
        singleton_correction: bool,

        /// Skip Stage 1 if BAM already exists in output/alignment/
        #[arg(long)]
        skip_alignment: bool,
    },

    /// Stage 1 only: extract UMI barcodes from FASTQ and align reads
    #[command(name = "fastq2bam")]
    Fastq2Bam {
        #[arg(long)]
        fastq1: PathBuf,
        #[arg(long)]
        fastq2: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(long)]
        bpattern: Option<String>,
        #[arg(long)]
        blist: Option<PathBuf>,
        #[arg(long, default_value = "bwa-mem2", value_parser = ["parabricks", "bwa-mem2", "minimap2", "bwa"])]
        aligner: String,
        #[arg(long)]
        aligner_path: Option<String>,
        #[arg(long)]
        reference: PathBuf,
        #[arg(long)]
        aligner_args: Option<String>,
        #[arg(long)]
        parabricks_docker: bool,
        #[arg(long, default_value = "nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1")]
        parabricks_image: String,
        #[arg(long, default_value = "samtools")]
        samtools: String,
    },

    /// Stage 2 only: generate consensus sequences from UMI-tagged BAM
    Consensus {
        /// Input BAM file (position-sorted, UMI in query name)
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(long)]
        bedfile: Option<PathBuf>,
        #[arg(long, default_value = "0.7")]
        cutoff: f64,
        #[arg(long, default_value = "30")]
        min_qual: u8,
        #[arg(long, default_value = "1")]
        min_family_size: usize,
        #[arg(long, default_value = "1000")]
        max_family_size: usize,
        #[arg(long)]
        singleton_correction: bool,
    },

    /// Generate a sample manifest from a directory of FASTQ files
    Manifest {
        /// Directory containing FASTQ files
        #[arg(short, long)]
        dir: PathBuf,

        /// Output manifest file (TSV)
        #[arg(short, long, default_value = "manifest.tsv")]
        output: PathBuf,

        /// Recursively scan subdirectories
        #[arg(long)]
        recursive: bool,
    },

    /// Batch process all samples in a manifest
    Batch {
        /// Manifest file (TSV: sample_name, fastq1, fastq2)
        /// Generate with `anneal manifest --dir /path/to/fastqs`
        #[arg(short, long)]
        manifest: PathBuf,

        /// Top-level output directory (each sample gets a subdirectory)
        #[arg(short, long)]
        output: PathBuf,

        /// Reference genome (FASTA, indexed)
        #[arg(long)]
        reference: PathBuf,

        /// Barcode pattern
        #[arg(long, default_value = "NNNSS")]
        bpattern: String,

        /// Aligner backend
        #[arg(long, default_value = "bwa-mem2", value_parser = ["parabricks", "bwa-mem2", "minimap2", "bwa"])]
        aligner: String,

        /// Path to aligner executable
        #[arg(long)]
        aligner_path: Option<String>,

        /// Extra arguments passed to the aligner
        #[arg(long)]
        aligner_args: Option<String>,

        /// For Parabricks: use Docker
        #[arg(long)]
        parabricks_docker: bool,

        /// For Parabricks: Docker image tag
        #[arg(long, default_value = "nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1")]
        parabricks_image: String,

        /// Path to samtools
        #[arg(long, default_value = "samtools")]
        samtools: String,

        /// BED file for target regions
        #[arg(long)]
        bedfile: Option<PathBuf>,

        /// Consensus cutoff fraction
        #[arg(long, default_value = "0.7")]
        cutoff: f64,

        /// Minimum base quality
        #[arg(long, default_value = "30")]
        min_qual: u8,

        /// Minimum reads per family
        #[arg(long, default_value = "1")]
        min_family_size: usize,

        /// Maximum family size
        #[arg(long, default_value = "1000")]
        max_family_size: usize,

        /// Enable singleton correction
        #[arg(long)]
        singleton_correction: bool,

        /// Skip Stage 1 for samples that already have a BAM
        #[arg(long)]
        skip_alignment: bool,
    },
}

// ---------------------------------------------------------------------------
// GPU initialization
// ---------------------------------------------------------------------------

fn init_gpu(no_gpu: bool, _gpu_device: usize) -> bool {
    if no_gpu {
        info!("GPU acceleration disabled by --no-gpu flag");
        return false;
    }
    #[cfg(feature = "gpu")]
    {
        match cuda::gpu::GpuContext::new(_gpu_device) {
            Ok(ctx) => {
                info!("GPU initialized: {}", ctx.device_name());
                true
            }
            Err(e) => {
                info!("GPU not available ({}), falling back to CPU", e);
                false
            }
        }
    }
    #[cfg(not(feature = "gpu"))]
    {
        info!("Built without GPU support, using CPU mode");
        false
    }
}

// ---------------------------------------------------------------------------
// Stage runners
// ---------------------------------------------------------------------------

fn run_stage1(
    fastq1: &PathBuf,
    fastq2: &PathBuf,
    output: &PathBuf,
    bpattern: Option<&str>,
    blist: Option<&PathBuf>,
    aligner: &str,
    aligner_path: Option<String>,
    reference: &PathBuf,
    aligner_args: Option<String>,
    parabricks_docker: bool,
    parabricks_image: &str,
    samtools: &str,
    threads: usize,
    gpu_device: usize,
) -> Result<PathBuf> {
    info!("=== Stage 1: Barcode extraction and alignment ===");
    info!("Aligner backend: {}", aligner);

    let aligner_config = barcode::extract::AlignerConfig {
        backend: aligner.parse()?,
        executable_path: aligner_path,
        reference: reference.clone(),
        extra_args: aligner_args,
        parabricks_docker,
        parabricks_image: parabricks_image.to_string(),
        samtools_path: samtools.to_string(),
        threads,
        gpu_device,
    };

    barcode::extract::run_fastq2bam(
        fastq1,
        fastq2,
        output,
        bpattern,
        blist.map(|p| p.as_path()),
        &aligner_config,
    )?;

    // Locate the output BAM
    let stem = fastq1
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("output");
    let bam_stem = stem
        .strip_suffix(".fastq")
        .or_else(|| stem.strip_suffix(".fq"))
        .unwrap_or(stem);
    let bam_path = output.join(format!("{}.bam", bam_stem));

    if bam_path.exists() {
        info!("Stage 1 output BAM: {}", bam_path.display());
        return Ok(bam_path);
    }

    // Fallback: find any BAM in output directory
    let mut bams: Vec<_> = std::fs::read_dir(output)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().map_or(false, |ext| ext == "bam"))
        .collect();
    bams.sort();
    bams.into_iter()
        .next()
        .context("No BAM file found in Stage 1 output directory")
}

fn run_stage2(
    input: &PathBuf,
    output: &PathBuf,
    bedfile: Option<&PathBuf>,
    cutoff: f64,
    min_qual: u8,
    min_family_size: usize,
    max_family_size: usize,
    singleton_correction: bool,
    gpu_available: bool,
    gpu_device: usize,
) -> Result<()> {
    info!("=== Stage 2: Consensus sequence generation ===");

    let config = consensus::config::ConsensusConfig {
        cutoff,
        min_qual,
        min_family_size,
        max_family_size,
        use_gpu: gpu_available,
        gpu_device,
        singleton_correction,
    };

    consensus::pipeline::run_consensus(input, output, bedfile.map(|p| p.as_path()), &config)
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Configure logging
    let log_level = match cli.verbose {
        0 => "info",
        1 => "debug",
        _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level)).init();

    // Set thread count
    let threads = cli.threads.unwrap_or_else(num_cpus::get);
    if let Some(t) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()?;
    }

    let gpu_available = init_gpu(cli.no_gpu, cli.gpu_device);
    let start = std::time::Instant::now();

    match cli.command {
        // =================================================================
        // Full pipeline
        // =================================================================
        Commands::Run {
            fastq1,
            fastq2,
            output,
            bpattern,
            blist,
            aligner,
            aligner_path,
            reference,
            aligner_args,
            parabricks_docker,
            parabricks_image,
            samtools,
            bedfile,
            cutoff,
            min_qual,
            min_family_size,
            max_family_size,
            singleton_correction,
            skip_alignment,
        } => {
            info!("=== Anneal v2: Full duplex consensus pipeline ===");

            let align_dir = output.join("alignment");
            let consensus_dir = output.join("consensus");
            std::fs::create_dir_all(&align_dir)?;
            std::fs::create_dir_all(&consensus_dir)?;

            // Stage 1: Alignment
            let bam_path = if skip_alignment {
                let mut bams: Vec<_> = std::fs::read_dir(&align_dir)?
                    .filter_map(|e| e.ok())
                    .map(|e| e.path())
                    .filter(|p| p.extension().map_or(false, |ext| ext == "bam"))
                    .collect();
                bams.sort();
                let bam = bams
                    .into_iter()
                    .next()
                    .context("--skip-alignment set but no BAM found in alignment/ directory")?;
                info!("Skipping alignment, using existing BAM: {}", bam.display());
                bam
            } else {
                let t1 = std::time::Instant::now();
                let bam = run_stage1(
                    &fastq1,
                    &fastq2,
                    &align_dir,
                    Some(bpattern.as_str()),
                    blist.as_ref(),
                    &aligner,
                    aligner_path,
                    &reference,
                    aligner_args,
                    parabricks_docker,
                    &parabricks_image,
                    &samtools,
                    threads,
                    cli.gpu_device,
                )?;
                info!("Stage 1 completed in {:.1}s", t1.elapsed().as_secs_f64());
                bam
            };

            // Stage 2: Consensus
            let t2 = std::time::Instant::now();
            run_stage2(
                &bam_path,
                &consensus_dir,
                bedfile.as_ref(),
                cutoff,
                min_qual,
                min_family_size,
                max_family_size,
                singleton_correction,
                gpu_available,
                cli.gpu_device,
            )?;
            info!("Stage 2 completed in {:.1}s", t2.elapsed().as_secs_f64());

            let elapsed = start.elapsed().as_secs_f64();
            info!("=== Pipeline complete: {:.1}s ({:.1} min) ===", elapsed, elapsed / 60.0);
            info!("Output:");
            info!("  Aligned BAM:  {}", bam_path.display());
            if singleton_correction {
                info!("  DCS BAM:      {}", consensus_dir.join("dcs.sc.sorted.bam").display());
                info!("  SSCS BAM:     {}", consensus_dir.join("sscs.sc.sorted.bam").display());
            } else {
                info!("  DCS BAM:      {}", consensus_dir.join("dcs.sorted.bam").display());
                info!("  SSCS BAM:     {}", consensus_dir.join("sscs.sorted.bam").display());
            }
            info!("  Statistics:   {}", consensus_dir.join("stats.txt").display());
        }

        // =================================================================
        // Stage 1 only
        // =================================================================
        Commands::Fastq2Bam {
            fastq1,
            fastq2,
            output,
            bpattern,
            blist,
            aligner,
            aligner_path,
            reference,
            aligner_args,
            parabricks_docker,
            parabricks_image,
            samtools,
        } => {
            run_stage1(
                &fastq1,
                &fastq2,
                &output,
                bpattern.as_deref(),
                blist.as_ref(),
                &aligner,
                aligner_path,
                &reference,
                aligner_args,
                parabricks_docker,
                &parabricks_image,
                &samtools,
                threads,
                cli.gpu_device,
            )?;
        }

        // =================================================================
        // Stage 2 only
        // =================================================================
        Commands::Consensus {
            input,
            output,
            bedfile,
            cutoff,
            min_qual,
            min_family_size,
            max_family_size,
            singleton_correction,
        } => {
            run_stage2(
                &input,
                &output,
                bedfile.as_ref(),
                cutoff,
                min_qual,
                min_family_size,
                max_family_size,
                singleton_correction,
                gpu_available,
                cli.gpu_device,
            )?;
        }

        // =================================================================
        // Manifest generation
        // =================================================================
        Commands::Manifest {
            dir,
            output,
            recursive,
        } => {
            info!("Scanning {} for paired FASTQ files...", dir.display());
            let samples = manifest::scan_directory(&dir, recursive)?;

            for (i, s) in samples.iter().enumerate() {
                info!(
                    "  [{}] {} | R1: {} | R2: {}",
                    i + 1,
                    s.sample_name,
                    s.fastq1.display(),
                    s.fastq2.display()
                );
            }

            manifest::write_manifest(&samples, &output)?;
            info!("Manifest with {} sample(s) written to {}", samples.len(), output.display());
        }

        // =================================================================
        // Batch processing
        // =================================================================
        Commands::Batch {
            manifest: manifest_path,
            output,
            reference,
            bpattern,
            aligner,
            aligner_path,
            aligner_args,
            parabricks_docker,
            parabricks_image,
            samtools,
            bedfile,
            cutoff,
            min_qual,
            min_family_size,
            max_family_size,
            singleton_correction,
            skip_alignment,
        } => {
            let samples = manifest::read_manifest(&manifest_path)?;
            let n_samples = samples.len();

            info!("=== Anneal v2: Batch processing {} sample(s) ===", n_samples);
            std::fs::create_dir_all(&output)?;

            // Copy manifest to output directory for reproducibility
            let manifest_copy = output.join("manifest.tsv");
            manifest::write_manifest(&samples, &manifest_copy)?;

            let mut completed = 0u32;
            let mut failed = Vec::new();

            for (i, sample) in samples.iter().enumerate() {
                let sample_start = std::time::Instant::now();
                info!(
                    "=== [{}/{}] Sample: {} ===",
                    i + 1,
                    n_samples,
                    sample.sample_name
                );

                // Validate FASTQ files exist
                if !sample.fastq1.exists() {
                    warn!("R1 not found: {} -- skipping", sample.fastq1.display());
                    failed.push(sample.sample_name.clone());
                    continue;
                }
                if !sample.fastq2.exists() {
                    warn!("R2 not found: {} -- skipping", sample.fastq2.display());
                    failed.push(sample.sample_name.clone());
                    continue;
                }

                // Per-sample output directory
                let sample_dir = output.join(&sample.sample_name);
                let align_dir = sample_dir.join("alignment");
                let consensus_dir = sample_dir.join("consensus");
                std::fs::create_dir_all(&align_dir)?;
                std::fs::create_dir_all(&consensus_dir)?;

                // Stage 1
                let bam_path = if skip_alignment {
                    let mut bams: Vec<_> = std::fs::read_dir(&align_dir)
                        .ok()
                        .map(|rd| {
                            rd.filter_map(|e| e.ok())
                                .map(|e| e.path())
                                .filter(|p| p.extension().map_or(false, |ext| ext == "bam"))
                                .collect()
                        })
                        .unwrap_or_default();
                    bams.sort();
                    match bams.into_iter().next() {
                        Some(bam) => {
                            info!("Skipping alignment, using: {}", bam.display());
                            bam
                        }
                        None => {
                            warn!(
                                "--skip-alignment but no BAM in {} -- running alignment",
                                align_dir.display()
                            );
                            match run_stage1(
                                &sample.fastq1,
                                &sample.fastq2,
                                &align_dir,
                                Some(bpattern.as_str()),
                                None,
                                &aligner,
                                aligner_path.clone(),
                                &reference,
                                aligner_args.clone(),
                                parabricks_docker,
                                &parabricks_image,
                                &samtools,
                                threads,
                                cli.gpu_device,
                            ) {
                                Ok(bam) => bam,
                                Err(e) => {
                                    warn!("Stage 1 failed for {}: {}", sample.sample_name, e);
                                    failed.push(sample.sample_name.clone());
                                    continue;
                                }
                            }
                        }
                    }
                } else {
                    match run_stage1(
                        &sample.fastq1,
                        &sample.fastq2,
                        &align_dir,
                        Some(bpattern.as_str()),
                        None,
                        &aligner,
                        aligner_path.clone(),
                        &reference,
                        aligner_args.clone(),
                        parabricks_docker,
                        &parabricks_image,
                        &samtools,
                        threads,
                        cli.gpu_device,
                    ) {
                        Ok(bam) => bam,
                        Err(e) => {
                            warn!("Stage 1 failed for {}: {}", sample.sample_name, e);
                            failed.push(sample.sample_name.clone());
                            continue;
                        }
                    }
                };

                // Stage 2
                match run_stage2(
                    &bam_path,
                    &consensus_dir,
                    bedfile.as_ref(),
                    cutoff,
                    min_qual,
                    min_family_size,
                    max_family_size,
                    singleton_correction,
                    gpu_available,
                    cli.gpu_device,
                ) {
                    Ok(()) => {
                        completed += 1;
                        info!(
                            "Sample {} complete in {:.1}s",
                            sample.sample_name,
                            sample_start.elapsed().as_secs_f64()
                        );
                    }
                    Err(e) => {
                        warn!("Stage 2 failed for {}: {}", sample.sample_name, e);
                        failed.push(sample.sample_name.clone());
                    }
                }
            }

            // Write batch summary
            let summary_path = output.join("batch_summary.txt");
            write_batch_summary(&summary_path, n_samples, completed, &failed, start.elapsed())?;

            let elapsed = start.elapsed().as_secs_f64();
            info!("=== Batch complete ===");
            info!("  Samples processed: {}/{}", completed, n_samples);
            if !failed.is_empty() {
                info!("  Failed: {}", failed.join(", "));
            }
            info!("  Total time: {:.1}s ({:.1} min)", elapsed, elapsed / 60.0);
            info!("  Summary: {}", summary_path.display());
        }
    }

    info!("Pipeline complete.");
    Ok(())
}

/// Write a batch processing summary.
fn write_batch_summary(
    path: &Path,
    total: usize,
    completed: u32,
    failed: &[String],
    elapsed: std::time::Duration,
) -> Result<()> {
    use std::io::Write;
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(f, "# Anneal Batch Summary")?;
    writeln!(f, "Date: {}", chrono_now())?;
    writeln!(f, "Total samples: {}", total)?;
    writeln!(f, "Completed: {}", completed)?;
    writeln!(f, "Failed: {}", failed.len())?;
    writeln!(
        f,
        "Total time: {:.1}s ({:.1} min)",
        elapsed.as_secs_f64(),
        elapsed.as_secs_f64() / 60.0
    )?;
    if !failed.is_empty() {
        writeln!(f, "\nFailed samples:")?;
        for s in failed {
            writeln!(f, "  {}", s)?;
        }
    }
    Ok(())
}

/// Simple timestamp without chrono dependency.
fn chrono_now() -> String {
    let output = std::process::Command::new("date")
        .arg("+%Y-%m-%d %H:%M:%S")
        .output();
    match output {
        Ok(o) => String::from_utf8_lossy(&o.stdout).trim().to_string(),
        Err(_) => "unknown".to_string(),
    }
}
