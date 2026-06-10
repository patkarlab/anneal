//! Duplex Consensus Sequence (DCS) generation.
//!
//! Pairs complementary SSCS reads and generates duplex consensus.
//! Pairing rule: swap UMI halves + flip strand.

use ahash::AHashMap;
use crate::consensus::config::ConsensusConfig;
use crate::consensus::cpu::call_duplex_consensus_cpu;
use crate::consensus::sscs::{SscsRead, ReadMetadata};
use crate::grouping::families::build_duplex_complement_tag;
use log::{debug, info};

/// Result of DCS generation.
pub struct DcsResult {
    pub dcs_reads: Vec<DcsRead>,
    pub sscs_singletons: Vec<SscsRead>,
}

/// A duplex consensus read.
#[derive(Debug)]
pub struct DcsRead {
    pub tag: String,
    pub consensus: crate::consensus::cpu::ConsensusResult,
    pub representative: ReadMetadata,
}

/// Generate DCS by pairing complementary SSCS reads.
pub fn generate_dcs(sscs_reads: Vec<SscsRead>, _config: &ConsensusConfig) -> DcsResult {
    debug!("Pairing {} SSCS reads for DCS", sscs_reads.len());

    // Index by tag
    let mut by_tag: AHashMap<String, SscsRead> = AHashMap::with_capacity(sscs_reads.len());
    for s in sscs_reads {
        by_tag.insert(s.tag.clone(), s);
    }

    let mut dcs_reads = Vec::new();
    let mut paired: ahash::AHashSet<String> = ahash::AHashSet::new();

    let tags: Vec<String> = by_tag.keys().cloned().collect();

    for tag in &tags {
        if paired.contains(tag.as_str()) {
            continue;
        }
        let comp = match build_duplex_complement_tag(tag) {
            Some(c) => c,
            None => continue,
        };
        if paired.contains(comp.as_str()) || !by_tag.contains_key(&comp) {
            continue;
        }

        // Mark both as paired
        paired.insert(tag.clone());
        paired.insert(comp.clone());

        // Call duplex consensus
        let a = by_tag.get(tag).unwrap();
        let b = by_tag.get(&comp).unwrap();
        let duplex = call_duplex_consensus_cpu(&a.consensus, &b.consensus);
        let dcs_tag = format!("{}:{}+{}:{}", tag, a.consensus.family_size, comp, b.consensus.family_size);

        dcs_reads.push(DcsRead {
            tag: dcs_tag,
            consensus: duplex,
            representative: a.representative.clone(),
        });
    }

    // Collect unpaired
    let singletons: Vec<SscsRead> = by_tag
        .into_iter()
        .filter(|(t, _)| !paired.contains(t.as_str()))
        .map(|(_, s)| s)
        .collect();

    info!("DCS: {} pairs, {} unpaired SSCS", dcs_reads.len(), singletons.len());

    DcsResult {
        dcs_reads,
        sscs_singletons: singletons,
    }
}
