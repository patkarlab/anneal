//! Single-Strand Consensus Sequence (SSCS) generation.
//!
//! Takes read families (grouped by UMI + genomic position) and produces
//! one consensus read per family. Families of size 1 (singletons) are
//! separated for optional singleton correction downstream.

use crate::consensus::config::ConsensusConfig;
use crate::consensus::cpu::{call_consensus_cpu, ConsensusResult};
use crate::grouping::families::{ReadFamily, ReadData};
use log::debug;
use rayon::prelude::*;

/// Result of SSCS generation for one family.
#[derive(Debug, Clone)]
pub struct SscsRead {
    /// The family tag (used later for duplex pairing)
    pub tag: String,
    /// Consensus sequence and quality
    pub consensus: ConsensusResult,
    /// Representative read metadata (flag, pos, cigar, etc.)
    pub representative: ReadMetadata,
}

/// Minimal metadata carried forward from the family to construct output BAM records.
#[derive(Debug, Clone)]
pub struct ReadMetadata {
    pub flag: u16,
    pub rname: String,
    pub ref_id: i32,
    pub pos: i64,
    pub mapq: u8,
    pub cigar_str: String,
    pub mate_rname: String,
    pub mate_ref_id: i32,
    pub mate_pos: i64,
    pub template_len: i64,
}

impl ReadMetadata {
    /// Build from a ReadData, using the most common flag/mapq from the family.
    pub fn from_reads(reads: &[ReadData]) -> Self {
        // Use the first read as the representative; in a full implementation
        // you would pick the most common flag and mapping quality.
        let r = &reads[0];
        Self {
            flag: most_common_flag(reads),
            rname: r.rname.clone(),
            ref_id: r.ref_id,
            pos: r.pos,
            mapq: most_common_mapq(reads),
            cigar_str: r.cigar_str.clone(),
            mate_rname: r.mate_rname.clone(),
            mate_ref_id: r.mate_ref_id,
            mate_pos: r.mate_pos,
            template_len: r.template_len,
        }
    }
}

/// Find the most common SAM flag among reads.
/// Prioritizes properly paired flags (99/83/147/163).
fn most_common_flag(reads: &[ReadData]) -> u16 {
    let mut counts = std::collections::HashMap::new();
    for r in reads {
        *counts.entry(r.flag).or_insert(0u32) += 1;
    }

    // Prefer canonical paired flags
    const PREFERRED: [u16; 4] = [99, 83, 147, 163];
    let max_count = counts.values().max().copied().unwrap_or(0);

    for &pf in &PREFERRED {
        if counts.get(&pf).copied().unwrap_or(0) == max_count {
            return pf;
        }
    }

    // Fall back to most frequent
    counts
        .into_iter()
        .max_by_key(|(_, c)| *c)
        .map(|(f, _)| f)
        .unwrap_or(0)
}

fn most_common_mapq(reads: &[ReadData]) -> u8 {
    let mut counts = std::collections::HashMap::new();
    for r in reads {
        *counts.entry(r.mapq).or_insert(0u32) += 1;
    }
    counts
        .into_iter()
        .max_by_key(|(_, c)| *c)
        .map(|(q, _)| q)
        .unwrap_or(0)
}

/// Generate SSCS from a collection of read families.
///
/// Uses rayon for parallel processing across families on CPU.
/// When GPU is enabled, large batches are dispatched to CUDA kernels.
pub fn generate_sscs(
    families: Vec<ReadFamily>,
    config: &ConsensusConfig,
) -> Vec<SscsRead> {
    debug!(
        "Generating SSCS for {} families (GPU: {})",
        families.len(),
        config.use_gpu
    );

    if config.use_gpu {
        // TODO: dispatch to GPU path for large batches
        // For now, fall through to CPU parallel path
        generate_sscs_cpu_parallel(families, config)
    } else {
        generate_sscs_cpu_parallel(families, config)
    }
}

/// CPU parallel SSCS generation using rayon.
fn generate_sscs_cpu_parallel(
    families: Vec<ReadFamily>,
    config: &ConsensusConfig,
) -> Vec<SscsRead> {
    families
        .into_par_iter()
        .map(|family| {
            let metadata = ReadMetadata::from_reads(&family.reads);
            let consensus = call_consensus_cpu(
                &family.reads,
                config.cutoff,
                config.min_qual,
            );

            SscsRead {
                tag: family.tag,
                consensus,
                representative: metadata,
            }
        })
        .collect()
}
