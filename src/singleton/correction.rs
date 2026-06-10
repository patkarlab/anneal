//! Singleton correction implementation.
//!
//! Two rescue strategies:
//!   1. Singleton + complementary SSCS (sscs.rescue)
//!      - Find the duplex complement tag for the singleton
//!      - If an SSCS exists for that tag, use it to correct the singleton
//!
//!   2. Singleton + complementary singleton (singleton.rescue)
//!      - If both strands have unpaired singletons with matching duplex
//!        complement tags, combine them

use crate::consensus::cpu::{call_duplex_consensus_cpu, ConsensusResult};
use crate::consensus::sscs::SscsRead;
use crate::grouping::families::{build_duplex_complement_tag, ReadFamily};
use ahash::AHashMap;
use log::{debug, info};

/// Results of singleton correction.
pub struct CorrectionResult {
    /// Singletons rescued by complementary SSCS
    pub sscs_rescued: Vec<SscsRead>,
    /// Singletons rescued by complementary singletons
    pub singleton_rescued: Vec<SscsRead>,
    /// Singletons that could not be rescued
    pub remaining: Vec<ReadFamily>,
}

/// Attempt to correct singletons using complementary strand information.
///
/// Strategy 1: For each singleton, check if an SSCS exists on the
///   complementary strand (swapped UMI, flipped orientation).
///   If found, perform duplex-like correction.
///
/// Strategy 2: For remaining singletons, check if a complementary
///   singleton exists. If found, combine them.
pub fn correct_singletons(
    singletons: Vec<ReadFamily>,
    sscs_reads: &[SscsRead],
) -> CorrectionResult {
    // Index SSCS reads by tag for fast lookup
    let sscs_index: AHashMap<&str, &SscsRead> = sscs_reads
        .iter()
        .map(|s| (s.tag.as_str(), s))
        .collect();

    let mut sscs_rescued = Vec::new();
    let mut remaining_singletons = Vec::new();

    // Strategy 1: Singleton + complementary SSCS
    for singleton in singletons {
        if let Some(complement_tag) = build_duplex_complement_tag(&singleton.tag) {
            if let Some(complement_sscs) = sscs_index.get(complement_tag.as_str()) {
                // Found complementary SSCS: correct the singleton using
                // the SSCS as the "truth" from the opposite strand
                debug!(
                    "Rescuing singleton {} with complementary SSCS {}",
                    singleton.tag, complement_tag
                );

                // Create a "consensus" from the single read
                let singleton_read = &singleton.reads[0];
                let singleton_consensus = ConsensusResult {
                    sequence: singleton_read.sequence.clone(),
                    qualities: singleton_read.qualities.clone(),
                    proportions: vec![1.0; singleton_read.sequence.len()],
                    family_size: 1,
                };

                // Duplex-correct using the SSCS
                let corrected = call_duplex_consensus_cpu(
                    &singleton_consensus,
                    &complement_sscs.consensus,
                );

                sscs_rescued.push(SscsRead {
                    tag: singleton.tag.clone(),
                    consensus: corrected,
                    representative: crate::consensus::sscs::ReadMetadata::from_reads(&singleton.reads),
                });

                continue;
            }
        }
        remaining_singletons.push(singleton);
    }

    // Strategy 2: Singleton + complementary singleton
    let mut singleton_index: AHashMap<String, ReadFamily> = AHashMap::new();
    let mut singleton_rescued = Vec::new();
    let mut final_remaining = Vec::new();
    let mut used_tags: ahash::AHashSet<String> = ahash::AHashSet::new();

    for singleton in remaining_singletons {
        singleton_index.insert(singleton.tag.clone(), singleton);
    }

    let tags: Vec<String> = singleton_index.keys().cloned().collect();

    for tag in &tags {
        if used_tags.contains(tag.as_str()) {
            continue;
        }

        if let Some(complement_tag) = build_duplex_complement_tag(tag) {
            if used_tags.contains(complement_tag.as_str()) {
                continue;
            }

            if singleton_index.contains_key(&complement_tag) {
                debug!(
                    "Rescuing singleton pair: {} + {}",
                    tag, complement_tag
                );
                used_tags.insert(tag.clone());
                used_tags.insert(complement_tag.clone());

                // Both singletons found: combine them
                if let (Some(s1), Some(s2)) = (
                    singleton_index.get(tag),
                    singleton_index.get(&complement_tag),
                ) {
                    let c1 = ConsensusResult {
                        sequence: s1.reads[0].sequence.clone(),
                        qualities: s1.reads[0].qualities.clone(),
                        proportions: vec![1.0; s1.reads[0].sequence.len()],
                        family_size: 1,
                    };
                    let c2 = ConsensusResult {
                        sequence: s2.reads[0].sequence.clone(),
                        qualities: s2.reads[0].qualities.clone(),
                        proportions: vec![1.0; s2.reads[0].sequence.len()],
                        family_size: 1,
                    };
                    let corrected = call_duplex_consensus_cpu(&c1, &c2);

                    singleton_rescued.push(SscsRead {
                        tag: tag.clone(),
                        consensus: corrected,
                        representative: crate::consensus::sscs::ReadMetadata::from_reads(&s1.reads),
                    });
                }
            }
        }
    }

    // Collect remaining uncorrected singletons
    for (tag, family) in singleton_index {
        if !used_tags.contains(&tag) {
            final_remaining.push(family);
        }
    }

    info!(
        "Singleton correction: {} by SSCS, {} by singletons, {} remaining",
        sscs_rescued.len(),
        singleton_rescued.len(),
        final_remaining.len(),
    );

    CorrectionResult {
        sscs_rescued,
        singleton_rescued,
        remaining: final_remaining,
    }
}
