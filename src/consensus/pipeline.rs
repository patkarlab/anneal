//! Consensus pipeline: full ConsensusCruncher-compatible output hierarchy.
//!
//! Output files (with singleton correction enabled):
//!   sscs.sorted.bam               - Standard SSCS from multi-read families
//!   singleton.sorted.bam          - Raw singleton reads (family size = 1)
//!   sscs.rescue.sorted.bam        - Singletons rescued by complementary SSCS (Strategy 1)
//!   singleton.rescue.sorted.bam   - Singletons rescued by complementary singleton (Strategy 2)
//!   sscs.sc.sorted.bam            - Expanded SSCS pool (RECOMMENDED for downstream)
//!   dcs.sc.sorted.bam             - DCS from expanded pool (RECOMMENDED for downstream)
//!   sscs.singleton.sorted.bam     - Unpaired SSCS after DCS formation
//!   rescue.remaining.sorted.bam   - Unrescuable singletons (no complement found)
//!   stats.txt                     - Comprehensive statistics
//!   family_sizes.tsv              - Family size distribution
//!
//! Without singleton correction:
//!   sscs.sorted.bam               - SSCS (RECOMMENDED)
//!   singleton.sorted.bam          - Raw singletons
//!   dcs.sorted.bam                - DCS (RECOMMENDED)
//!   sscs.singleton.sorted.bam     - Unpaired SSCS
//!   stats.txt
//!   family_sizes.tsv

use crate::consensus::config::ConsensusConfig;
use crate::consensus::{dcs, sscs};
use crate::consensus::sscs::SscsRead;
use crate::grouping::families::{self, ReadData};

use anyhow::{bail, Context, Result};
use log::info;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::{Command, Stdio};
use std::time::Instant;

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

pub fn run_consensus(
    input_bam: &Path,
    output_dir: &Path,
    bedfile: Option<&Path>,
    config: &ConsensusConfig,
) -> Result<()> {
    let start = Instant::now();
    std::fs::create_dir_all(output_dir)?;

    info!("Anneal 0.1.0 Consensus Pipeline");
    info!("Input BAM:  {}", input_bam.display());
    info!("Output dir: {}", output_dir.display());
    info!("Config: cutoff={}, min_qual={}, min_family={}, max_family={}, SC={}",
        config.cutoff, config.min_qual, config.min_family_size,
        config.max_family_size, config.singleton_correction);

    // Parse BAM header
    let header_text = get_sam_header(input_bam)?;

    // Determine regions
    let regions = if let Some(bed) = bedfile {
        info!("Using BED regions from: {}", bed.display());
        parse_bed(bed)?
    } else {
        let chroms = get_chromosomes_with_reads(input_bam)?;
        info!("Processing all {} chromosomes with mapped reads", chroms.len());
        chroms
            .into_iter()
            .map(|chrom| Region { chrom, start: None, end: None })
            .collect()
    };
    info!("Will process {} region(s)", regions.len());

    // -----------------------------------------------------------------------
    // Open SAM writers
    // -----------------------------------------------------------------------
    let sc = config.singleton_correction;

    // Always produced
    let sscs_sam = output_dir.join("sscs.sam");
    let singleton_sam = output_dir.join("singleton.sam");
    let unpaired_sam = output_dir.join("sscs.singleton.sam");
    let dcs_sam = output_dir.join(if sc { "dcs.sc.sam" } else { "dcs.sam" });

    let mut sscs_w = open_sam_writer(&sscs_sam, &header_text)?;
    let mut singleton_w = open_sam_writer(&singleton_sam, &header_text)?;
    let mut unpaired_w = open_sam_writer(&unpaired_sam, &header_text)?;
    let mut dcs_w = open_sam_writer(&dcs_sam, &header_text)?;

    // SC-specific writers
    let sscs_rescue_sam = output_dir.join("sscs.rescue.sam");
    let singleton_rescue_sam = output_dir.join("singleton.rescue.sam");
    let sscs_sc_sam = output_dir.join("sscs.sc.sam");
    let rescue_remaining_sam = output_dir.join("rescue.remaining.sam");

    let mut sscs_rescue_w = if sc { Some(open_sam_writer(&sscs_rescue_sam, &header_text)?) } else { None };
    let mut singleton_rescue_w = if sc { Some(open_sam_writer(&singleton_rescue_sam, &header_text)?) } else { None };
    let mut sscs_sc_w = if sc { Some(open_sam_writer(&sscs_sc_sam, &header_text)?) } else { None };
    let mut rescue_remaining_w = if sc { Some(open_sam_writer(&rescue_remaining_sam, &header_text)?) } else { None };

    // -----------------------------------------------------------------------
    // Statistics
    // -----------------------------------------------------------------------
    let mut stats = Stats::default();
    let mut family_sizes: BTreeMap<usize, u64> = BTreeMap::new();

    // -----------------------------------------------------------------------
    // Process each region
    // -----------------------------------------------------------------------
    for (i, region) in regions.iter().enumerate() {
        let t = Instant::now();
        let region_str = match (region.start, region.end) {
            (Some(s), Some(e)) => format!("{}:{}-{}", region.chrom, s, e),
            _ => region.chrom.clone(),
        };
        info!("[{}/{}] {}", i + 1, regions.len(), region_str);

        // Read aligned reads
        let reads = read_region(input_bam, region)?;
        let nr = reads.len() as u64;
        stats.total_reads += nr;
        if reads.is_empty() {
            continue;
        }

        // Group into families
        let families = families::group_reads_into_families(
            reads.into_iter(),
            config.max_family_size,
        );

        // Track family size distribution
        let nf = families.len();
        stats.total_families += nf as u64;
        for (_, fam) in &families {
            *family_sizes.entry(fam.family_size()).or_insert(0) += 1;
        }

        // Partition: multi-read families vs singletons
        let (multi, singles) = families::partition_families(families, config.min_family_size);
        stats.sscs_families += multi.len() as u64;

        // --- Write raw singletons ---
        for fam in &singles {
            write_read_sam(&mut singleton_w, &fam.reads[0])?;
        }
        stats.singletons += singles.len() as u64;

        // --- Generate standard SSCS ---
        let standard_sscs = sscs::generate_sscs(multi, config);
        stats.sscs_reads += standard_sscs.len() as u64;

        // Write standard SSCS to sscs.sam
        for s in &standard_sscs {
            write_sscs_sam(&mut sscs_w, s)?;
        }

        // If SC enabled, also write standard SSCS into the expanded pool file
        if let Some(ref mut w) = sscs_sc_w {
            for s in &standard_sscs {
                write_sscs_sam(w, s)?;
            }
        }

        // --- Build consensus pool for DCS ---
        let mut pool = standard_sscs;

        if sc && !singles.is_empty() {
            let cr = crate::singleton::correction::correct_singletons(singles, &pool);

            // Strategy 1: rescued by complementary SSCS
            stats.sscs_rescued += cr.sscs_rescued.len() as u64;
            if let Some(ref mut w) = sscs_rescue_w {
                for s in &cr.sscs_rescued {
                    write_sscs_sam(w, s)?;
                }
            }
            if let Some(ref mut w) = sscs_sc_w {
                for s in &cr.sscs_rescued {
                    write_sscs_sam(w, s)?;
                }
            }

            // Strategy 2: rescued by complementary singleton
            stats.singleton_rescued += cr.singleton_rescued.len() as u64;
            if let Some(ref mut w) = singleton_rescue_w {
                for s in &cr.singleton_rescued {
                    write_sscs_sam(w, s)?;
                }
            }
            if let Some(ref mut w) = sscs_sc_w {
                for s in &cr.singleton_rescued {
                    write_sscs_sam(w, s)?;
                }
            }

            // Remaining: unrescuable singletons
            stats.rescue_remaining += cr.remaining.len() as u64;
            if let Some(ref mut w) = rescue_remaining_w {
                for fam in &cr.remaining {
                    write_read_sam(w, &fam.reads[0])?;
                }
            }

            pool.extend(cr.sscs_rescued);
            pool.extend(cr.singleton_rescued);
        } else if !sc {
            // SC disabled: singletons are just counted, already written
        }

        stats.expanded_pool += pool.len() as u64;

        // --- Generate DCS ---
        let dcs_result = dcs::generate_dcs(pool, config);
        stats.dcs_reads += dcs_result.dcs_reads.len() as u64;
        stats.unpaired_sscs += dcs_result.sscs_singletons.len() as u64;

        // Write DCS
        for d in &dcs_result.dcs_reads {
            write_dcs_sam(&mut dcs_w, d)?;
        }

        // Write unpaired SSCS (those that couldn't form DCS)
        for s in &dcs_result.sscs_singletons {
            write_sscs_sam(&mut unpaired_w, s)?;
        }

        info!(
            "  {:.1}s | {} reads | {} families | {} SSCS | {} DCS",
            t.elapsed().as_secs_f64(), nr, nf,
            dcs_result.sscs_singletons.len() + dcs_result.dcs_reads.len(),
            dcs_result.dcs_reads.len(),
        );
    }

    // -----------------------------------------------------------------------
    // Flush and close all SAM writers
    // -----------------------------------------------------------------------
    sscs_w.flush()?;
    singleton_w.flush()?;
    unpaired_w.flush()?;
    dcs_w.flush()?;
    drop(sscs_w);
    drop(singleton_w);
    drop(unpaired_w);
    drop(dcs_w);

    if let Some(mut w) = sscs_rescue_w.take() { w.flush()?; }
    if let Some(mut w) = singleton_rescue_w.take() { w.flush()?; }
    if let Some(mut w) = sscs_sc_w.take() { w.flush()?; }
    if let Some(mut w) = rescue_remaining_w.take() { w.flush()?; }

    // -----------------------------------------------------------------------
    // Convert SAM -> sorted BAM + index
    // -----------------------------------------------------------------------
    info!("Converting SAM files to sorted BAM...");

    // Always produced
    convert_and_cleanup(&sscs_sam, output_dir, "sscs.sorted.bam")?;
    convert_and_cleanup(&singleton_sam, output_dir, "singleton.sorted.bam")?;
    convert_and_cleanup(&unpaired_sam, output_dir, "sscs.singleton.sorted.bam")?;

    let dcs_bam_name = if sc { "dcs.sc.sorted.bam" } else { "dcs.sorted.bam" };
    convert_and_cleanup(&dcs_sam, output_dir, dcs_bam_name)?;

    // SC-specific
    if sc {
        convert_and_cleanup(&sscs_rescue_sam, output_dir, "sscs.rescue.sorted.bam")?;
        convert_and_cleanup(&singleton_rescue_sam, output_dir, "singleton.rescue.sorted.bam")?;
        convert_and_cleanup(&sscs_sc_sam, output_dir, "sscs.sc.sorted.bam")?;
        convert_and_cleanup(&rescue_remaining_sam, output_dir, "rescue.remaining.sorted.bam")?;
    }

    // -----------------------------------------------------------------------
    // Write statistics
    // -----------------------------------------------------------------------
    write_stats(&output_dir.join("stats.txt"), &stats, config)?;
    write_family_sizes(&output_dir.join("family_sizes.tsv"), &family_sizes)?;

    // -----------------------------------------------------------------------
    // Summary
    // -----------------------------------------------------------------------
    let elapsed = start.elapsed();
    info!("=== Consensus Pipeline Summary ===");
    info!("Total reads:           {}", stats.total_reads);
    info!("Total families:        {}", stats.total_families);
    info!("  Multi-read (SSCS):   {}", stats.sscs_families);
    info!("  Singletons:          {}", stats.singletons);
    info!("SSCS generated:        {}", stats.sscs_reads);
    if sc {
        info!("Singleton Correction:");
        info!("  Rescued by SSCS:     {} (Strategy 1)", stats.sscs_rescued);
        info!("  Rescued by singleton:{} (Strategy 2)", stats.singleton_rescued);
        info!("  Unrescuable:         {}", stats.rescue_remaining);
        info!("Expanded SSCS pool:    {}", stats.expanded_pool);
    }
    info!("DCS reads:             {}", stats.dcs_reads);
    info!("Unpaired SSCS:         {}", stats.unpaired_sscs);
    info!("");

    // Efficiency metrics
    if stats.total_reads > 0 {
        let sscs_eff = if sc {
            (stats.expanded_pool as f64 / stats.total_families as f64) * 100.0
        } else {
            (stats.sscs_reads as f64 / stats.total_families as f64) * 100.0
        };
        let dcs_eff = (stats.dcs_reads as f64 * 2.0 / stats.total_families as f64) * 100.0;
        info!("SSCS efficiency:       {:.1}%", sscs_eff);
        info!("DCS efficiency:        {:.1}%", dcs_eff);
    }

    info!("Time: {:.1}s ({:.1} min)", elapsed.as_secs_f64(), elapsed.as_secs_f64() / 60.0);
    info!("");
    if sc {
        info!("Recommended outputs:");
        info!("  SSCS: {}", output_dir.join("sscs.sc.sorted.bam").display());
        info!("  DCS:  {}", output_dir.join("dcs.sc.sorted.bam").display());
    } else {
        info!("Recommended outputs:");
        info!("  SSCS: {}", output_dir.join("sscs.sorted.bam").display());
        info!("  DCS:  {}", output_dir.join("dcs.sorted.bam").display());
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// SAM writer helpers
// ---------------------------------------------------------------------------

fn open_sam_writer(path: &Path, header: &str) -> Result<BufWriter<File>> {
    let mut w = BufWriter::new(File::create(path)?);
    write!(w, "{}", header)?;
    Ok(w)
}

fn convert_and_cleanup(sam: &Path, output_dir: &Path, bam_name: &str) -> Result<()> {
    let bam = output_dir.join(bam_name);
    sam_to_sorted_bam(sam, &bam)?;
    let _ = std::fs::remove_file(sam);
    Ok(())
}

// ---------------------------------------------------------------------------
// SAM output: SSCS reads
// ---------------------------------------------------------------------------

fn write_sscs_sam(w: &mut impl Write, sscs: &SscsRead) -> Result<()> {
    let seq: String = sscs.consensus.sequence.iter().map(|&b| b as char).collect();
    let qual: String = sscs.consensus.qualities
        .iter()
        .map(|&q| (q.min(93) + 33) as char)
        .collect();

    let rname = &sscs.representative.rname;
    let mrnm = &sscs.representative.mate_rname;
    let mate_field = if mrnm == rname { "=" } else { mrnm.as_str() };

    writeln!(
        w,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tXV:Z:SSCS\tXW:i:{}",
        sscs.tag,
        sscs.representative.flag,
        rname,
        sscs.representative.pos + 1,
        sscs.representative.mapq,
        sscs.representative.cigar_str,
        mate_field,
        sscs.representative.mate_pos + 1,
        sscs.representative.template_len,
        seq,
        qual,
        sscs.consensus.family_size,
    )?;
    Ok(())
}

// ---------------------------------------------------------------------------
// SAM output: DCS reads
// ---------------------------------------------------------------------------

fn write_dcs_sam(w: &mut impl Write, dcs: &dcs::DcsRead) -> Result<()> {
    let seq: String = dcs.consensus.sequence.iter().map(|&b| b as char).collect();
    let qual: String = dcs.consensus.qualities
        .iter()
        .map(|&q| (q.min(93) + 33) as char)
        .collect();

    let rname = &dcs.representative.rname;
    let mrnm = &dcs.representative.mate_rname;
    let mate_field = if mrnm == rname { "=" } else { mrnm.as_str() };

    writeln!(
        w,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tXV:Z:DCS\tXW:i:{}",
        dcs.tag,
        dcs.representative.flag,
        rname,
        dcs.representative.pos + 1,
        dcs.representative.mapq,
        dcs.representative.cigar_str,
        mate_field,
        dcs.representative.mate_pos + 1,
        dcs.representative.template_len,
        seq,
        qual,
        dcs.consensus.family_size,
    )?;
    Ok(())
}

// ---------------------------------------------------------------------------
// SAM output: raw reads (singletons, unrescuable)
// ---------------------------------------------------------------------------

fn write_read_sam(w: &mut impl Write, read: &ReadData) -> Result<()> {
    let seq: String = read.sequence.iter().map(|&b| b as char).collect();
    let qual: String = read.qualities
        .iter()
        .map(|&q| (q.min(93) + 33) as char)
        .collect();

    let mate_field = if read.mate_rname == read.rname { "=" } else { read.mate_rname.as_str() };

    writeln!(
        w,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        read.qname,
        read.flag,
        read.rname,
        read.pos + 1,
        read.mapq,
        read.cigar_str,
        mate_field,
        read.mate_pos + 1,
        read.template_len,
        seq,
        qual,
    )?;
    Ok(())
}

// ---------------------------------------------------------------------------
// BAM reading via samtools
// ---------------------------------------------------------------------------

fn get_sam_header(bam: &Path) -> Result<String> {
    let out = Command::new("samtools")
        .args(["view", "-H"])
        .arg(bam)
        .output()
        .context("Failed to run samtools")?;
    if !out.status.success() {
        bail!("samtools view -H failed");
    }
    Ok(String::from_utf8(out.stdout)?)
}

fn parse_sq_names(header: &str) -> Vec<String> {
    header
        .lines()
        .filter(|l| l.starts_with("@SQ"))
        .filter_map(|l| {
            l.split('\t')
                .find(|f| f.starts_with("SN:"))
                .map(|f| f[3..].to_string())
        })
        .collect()
}

fn get_chromosomes_with_reads(bam: &Path) -> Result<Vec<String>> {
    let out = Command::new("samtools")
        .args(["idxstats"])
        .arg(bam)
        .output()
        .context("samtools idxstats failed")?;
    if !out.status.success() {
        bail!("samtools idxstats failed");
    }
    let mut chroms = Vec::new();
    for line in String::from_utf8(out.stdout)?.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() >= 3 {
            let mapped: u64 = f[2].parse().unwrap_or(0);
            if mapped > 0 && f[0] != "*" {
                chroms.push(f[0].to_string());
            }
        }
    }
    Ok(chroms)
}

#[derive(Debug)]
struct Region {
    chrom: String,
    start: Option<u64>,
    end: Option<u64>,
}

fn parse_bed(bed: &Path) -> Result<Vec<Region>> {
    let mut regions = Vec::new();
    for line in BufReader::new(File::open(bed)?).lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 { continue; }
        regions.push(Region {
            chrom: f[0].to_string(),
            start: Some(f[1].parse()?),
            end: Some(f[2].parse()?),
        });
    }
    Ok(regions)
}

fn read_region(bam: &Path, region: &Region) -> Result<Vec<ReadData>> {
    let mut cmd = Command::new("samtools");
    cmd.args(["view", "-F", "2820"]);
    cmd.arg(bam);

    match (region.start, region.end) {
        (Some(s), Some(e)) => {
            cmd.arg(format!("{}:{}-{}", region.chrom, s + 1, e));
        }
        _ => {
            cmd.arg(&region.chrom);
        }
    }

    cmd.stdout(Stdio::piped()).stderr(Stdio::null());

    let child = cmd.spawn().context("Cannot start samtools view")?;
    let reader = BufReader::with_capacity(8 * 1024 * 1024, child.stdout.unwrap());
    let delim = '|';
    let mut reads = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('@') { continue; }

        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 11 { continue; }

        let qname = f[0];
        let flag: u16 = match f[1].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        if flag & 0x1 == 0 { continue; }
        if !qname.contains(delim) { continue; }

        let rname = f[2];
        if rname == "*" { continue; }

        let pos: i64 = f[3].parse::<i64>().unwrap_or(0) - 1;
        let mapq: u8 = f[4].parse().unwrap_or(0);
        let cigar = f[5];
        let mrnm = if f[6] == "=" { rname } else { f[6] };
        let mpos: i64 = f[7].parse::<i64>().unwrap_or(0) - 1;
        let tlen: i64 = f[8].parse().unwrap_or(0);
        let seq = f[9];
        let qual_ascii = f[10];

        let tag = families::build_family_tag(qname, rname, pos, mrnm, mpos, cigar, flag, delim);
        let qualities: Vec<u8> = qual_ascii.bytes().map(|b| b.saturating_sub(33)).collect();

        reads.push(ReadData {
            qname: qname.to_string(),
            tag,
            sequence: seq.as_bytes().to_vec(),
            qualities,
            flag,
            rname: rname.to_string(),
            ref_id: 0,
            pos,
            mapq,
            cigar_str: cigar.to_string(),
            mate_rname: mrnm.to_string(),
            mate_ref_id: 0,
            mate_pos: mpos,
            template_len: tlen,
        });
    }

    Ok(reads)
}

// ---------------------------------------------------------------------------
// SAM -> sorted BAM
// ---------------------------------------------------------------------------

fn sam_to_sorted_bam(sam: &Path, bam: &Path) -> Result<()> {
    let view = Command::new("samtools")
        .args(["view", "-bS"])
        .arg(sam)
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .context("samtools view failed")?;

    let sort = Command::new("samtools")
        .args(["sort", "-o"])
        .arg(bam)
        .stdin(view.stdout.unwrap())
        .status()
        .context("samtools sort failed")?;

    if !sort.success() {
        bail!("samtools sort failed");
    }

    let _ = Command::new("samtools")
        .args(["index"])
        .arg(bam)
        .status();

    Ok(())
}

// ---------------------------------------------------------------------------
// Statistics
// ---------------------------------------------------------------------------

#[derive(Debug, Default)]
struct Stats {
    total_reads: u64,
    total_families: u64,
    sscs_families: u64,
    sscs_reads: u64,
    singletons: u64,
    sscs_rescued: u64,
    singleton_rescued: u64,
    rescue_remaining: u64,
    expanded_pool: u64,
    dcs_reads: u64,
    unpaired_sscs: u64,
}

fn write_stats(path: &Path, s: &Stats, config: &ConsensusConfig) -> Result<()> {
    let mut f = File::create(path)?;
    writeln!(f, "# Anneal 0.1.0 Consensus Statistics")?;
    writeln!(f, "#")?;
    writeln!(f, "# Parameters:")?;
    writeln!(f, "#   cutoff:                {}", config.cutoff)?;
    writeln!(f, "#   min_qual:              {}", config.min_qual)?;
    writeln!(f, "#   min_family_size:       {}", config.min_family_size)?;
    writeln!(f, "#   max_family_size:       {}", config.max_family_size)?;
    writeln!(f, "#   singleton_correction:  {}", config.singleton_correction)?;
    writeln!(f, "#")?;
    writeln!(f, "Total reads:               {}", s.total_reads)?;
    writeln!(f, "Total families:            {}", s.total_families)?;
    writeln!(f, "  Multi-read families:     {}", s.sscs_families)?;
    writeln!(f, "  Singletons:              {}", s.singletons)?;
    writeln!(f, "SSCS reads generated:      {}", s.sscs_reads)?;

    if config.singleton_correction {
        writeln!(f, "")?;
        writeln!(f, "Singleton Correction:")?;
        writeln!(f, "  Rescued by SSCS:         {} (Strategy 1)", s.sscs_rescued)?;
        writeln!(f, "  Rescued by singleton:    {} (Strategy 2)", s.singleton_rescued)?;
        writeln!(f, "  Total rescued:           {}", s.sscs_rescued + s.singleton_rescued)?;
        writeln!(f, "  Unrescuable remaining:   {}", s.rescue_remaining)?;
        writeln!(f, "")?;
        writeln!(f, "Expanded SSCS pool:        {}", s.expanded_pool)?;
    }

    writeln!(f, "")?;
    writeln!(f, "DCS reads:                 {}", s.dcs_reads)?;
    writeln!(f, "Unpaired SSCS:             {}", s.unpaired_sscs)?;

    // Efficiency metrics
    if s.total_families > 0 {
        writeln!(f, "")?;
        writeln!(f, "Efficiency:")?;

        let sscs_input = if config.singleton_correction {
            s.expanded_pool
        } else {
            s.sscs_reads
        };
        let sscs_eff = (sscs_input as f64 / s.total_families as f64) * 100.0;
        let dcs_eff = (s.dcs_reads as f64 * 2.0 / s.total_families as f64) * 100.0;
        let singleton_pct = (s.singletons as f64 / s.total_families as f64) * 100.0;

        writeln!(f, "  Singleton rate:          {:.1}%", singleton_pct)?;
        writeln!(f, "  SSCS efficiency:         {:.1}%", sscs_eff)?;
        writeln!(f, "  DCS efficiency:          {:.1}%", dcs_eff)?;

        if config.singleton_correction && s.sscs_reads > 0 {
            let dcs_recovery = (s.dcs_reads as f64) / (sscs_input as f64 / 2.0) * 100.0;
            writeln!(f, "  DCS recovery:            {:.1}%", dcs_recovery)?;
        }
    }

    writeln!(f, "")?;
    if config.singleton_correction {
        writeln!(f, "Recommended outputs:")?;
        writeln!(f, "  SSCS: sscs.sc.sorted.bam")?;
        writeln!(f, "  DCS:  dcs.sc.sorted.bam")?;
    } else {
        writeln!(f, "Recommended outputs:")?;
        writeln!(f, "  SSCS: sscs.sorted.bam")?;
        writeln!(f, "  DCS:  dcs.sorted.bam")?;
    }

    Ok(())
}

fn write_family_sizes(path: &Path, sizes: &BTreeMap<usize, u64>) -> Result<()> {
    let mut f = File::create(path)?;
    writeln!(f, "family_size\tcount")?;
    for (size, count) in sizes {
        writeln!(f, "{}\t{}", size, count)?;
    }
    info!("Family size distribution written to: {}", path.display());
    Ok(())
}
