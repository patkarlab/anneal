//! CPU-based consensus sequence generation.
//!
//! Implements the majority-rules consensus algorithm:
//!   - At each position, count bases with quality >= min_qual
//!   - If the most frequent base exceeds the cutoff fraction, call it
//!   - Otherwise, call 'N'
//!   - Consensus quality = sum of qualities from agreeing bases (capped at Q60)

use crate::grouping::families::ReadData;

/// Base encoding for efficient counting.
/// A=0, C=1, G=2, T=3, N=4
#[inline]
pub fn encode_base(b: u8) -> usize {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4, // N or any ambiguity
    }
}

/// Decode base index back to ASCII.
#[inline]
pub fn decode_base(idx: usize) -> u8 {
    match idx {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'N',
    }
}

/// Result of consensus calling for a single family.
#[derive(Debug, Clone)]
pub struct ConsensusResult {
    /// Consensus sequence
    pub sequence: Vec<u8>,
    /// Consensus quality scores (Phred)
    pub qualities: Vec<u8>,
    /// Proportion of agreeing bases at each position
    pub proportions: Vec<f64>,
    /// Number of reads used
    pub family_size: usize,
}

/// Generate a consensus sequence from a family of reads (CPU path).
///
/// Algorithm (majority-rules consensus per Lam et al.):
///   1. For each position, count bases with quality >= min_qual
///   2. Most frequent base must represent >= cutoff fraction of total
///   3. Consensus quality = sum of Phred scores from agreeing bases
///   4. Quality capped at 60 to avoid overflow in downstream tools
pub fn call_consensus_cpu(
    reads: &[ReadData],
    cutoff: f64,
    min_qual: u8,
) -> ConsensusResult {
    if reads.is_empty() {
        return ConsensusResult {
            sequence: Vec::new(),
            qualities: Vec::new(),
            proportions: Vec::new(),
            family_size: 0,
        };
    }

    // For a single read, return it directly
    if reads.len() == 1 {
        return ConsensusResult {
            sequence: reads[0].sequence.clone(),
            qualities: reads[0].qualities.clone(),
            proportions: vec![1.0; reads[0].sequence.len()],
            family_size: 1,
        };
    }

    // Determine consensus read length (use the most common length)
    let read_len = reads
        .iter()
        .map(|r| r.sequence.len())
        .max()
        .unwrap_or(0);

    let mut consensus_seq = Vec::with_capacity(read_len);
    let mut consensus_qual = Vec::with_capacity(read_len);
    let mut proportions = Vec::with_capacity(read_len);

    for pos in 0..read_len {
        // Count bases at this position (only those with quality >= min_qual)
        let mut counts = [0u32; 5]; // A, C, G, T, N
        let mut qual_sums = [0u64; 5];
        let mut total_valid = 0u32;

        for read in reads {
            if pos < read.sequence.len() && pos < read.qualities.len() {
                let qual = read.qualities[pos];
                if qual >= min_qual {
                    let base_idx = encode_base(read.sequence[pos]);
                    counts[base_idx] += 1;
                    qual_sums[base_idx] += qual as u64;
                    total_valid += 1;
                }
            }
        }

        if total_valid == 0 {
            consensus_seq.push(b'N');
            consensus_qual.push(0);
            proportions.push(0.0);
            continue;
        }

        // Find the most frequent base (excluding N)
        let (best_idx, best_count) = counts[..4]
            .iter()
            .enumerate()
            .max_by_key(|(_, &c)| c)
            .unwrap();

        let proportion = *best_count as f64 / total_valid as f64;

        if proportion >= cutoff {
            consensus_seq.push(decode_base(best_idx));

            // Consensus quality: sum of quality scores from agreeing bases,
            // capped at 60 to prevent issues with downstream tools.
            let qual_sum = qual_sums[best_idx].min(60) as u8;
            consensus_qual.push(qual_sum);
            proportions.push(proportion);
        } else {
            // No consensus at this position
            consensus_seq.push(b'N');
            consensus_qual.push(0);
            proportions.push(proportion);
        }
    }

    ConsensusResult {
        sequence: consensus_seq,
        qualities: consensus_qual,
        proportions,
        family_size: reads.len(),
    }
}

/// Generate a duplex consensus from two SSCS reads (complementary strands).
///
/// At each position:
///   - If both strands agree, use that base with summed quality
///   - If they disagree, call N with quality 0
///   - If one strand has N, use the other strand's base
pub fn call_duplex_consensus_cpu(
    sscs_pos: &ConsensusResult,
    sscs_neg: &ConsensusResult,
) -> ConsensusResult {
    let len = sscs_pos.sequence.len().max(sscs_neg.sequence.len());
    let mut sequence = Vec::with_capacity(len);
    let mut qualities = Vec::with_capacity(len);
    let mut proportions = Vec::with_capacity(len);

    for i in 0..len {
        let base_pos = sscs_pos.sequence.get(i).copied().unwrap_or(b'N');
        let base_neg = sscs_neg.sequence.get(i).copied().unwrap_or(b'N');
        let qual_pos = sscs_pos.qualities.get(i).copied().unwrap_or(0);
        let qual_neg = sscs_neg.qualities.get(i).copied().unwrap_or(0);

        match (base_pos, base_neg) {
            (b'N', b'N') => {
                sequence.push(b'N');
                qualities.push(0);
                proportions.push(0.0);
            }
            (b'N', b) => {
                sequence.push(b);
                qualities.push(qual_neg);
                proportions.push(0.5);
            }
            (b, b'N') => {
                sequence.push(b);
                qualities.push(qual_pos);
                proportions.push(0.5);
            }
            (a, b) if a == b => {
                sequence.push(a);
                // Sum qualities, cap at 60
                let combined_qual = (qual_pos as u16 + qual_neg as u16).min(60) as u8;
                qualities.push(combined_qual);
                proportions.push(1.0);
            }
            (_, _) => {
                // Disagreement between strands
                sequence.push(b'N');
                qualities.push(0);
                proportions.push(0.0);
            }
        }
    }

    ConsensusResult {
        sequence,
        qualities,
        proportions,
        family_size: sscs_pos.family_size + sscs_neg.family_size,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_read(seq: &[u8], qual: &[u8]) -> ReadData {
        ReadData {
            qname: "test".into(),
            tag: "test_tag".into(),
            sequence: seq.to_vec(),
            qualities: qual.to_vec(),
            flag: 99,
            rname: "chr1".into(),
            ref_id: 0,
            pos: 100,
            mapq: 60,
            cigar_str: "10M".into(),
            mate_rname: "chr1".into(),
            mate_ref_id: 0,
            mate_pos: 200,
            template_len: 110,
        }
    }

    #[test]
    fn test_consensus_unanimous() {
        let reads = vec![
            make_read(b"ACGTACGT", &[40, 40, 40, 40, 40, 40, 40, 40]),
            make_read(b"ACGTACGT", &[40, 40, 40, 40, 40, 40, 40, 40]),
            make_read(b"ACGTACGT", &[40, 40, 40, 40, 40, 40, 40, 40]),
        ];
        let result = call_consensus_cpu(&reads, 0.7, 30);
        assert_eq!(result.sequence, b"ACGTACGT");
        assert_eq!(result.family_size, 3);
    }

    #[test]
    fn test_consensus_majority() {
        let reads = vec![
            make_read(b"ACGT", &[40, 40, 40, 40]),
            make_read(b"ACGT", &[40, 40, 40, 40]),
            make_read(b"TCGT", &[40, 40, 40, 40]), // first base differs
        ];
        let result = call_consensus_cpu(&reads, 0.7, 30);
        // Position 0: A=2/3=0.67 < 0.7, so N
        assert_eq!(result.sequence[0], b'N');
        // Other positions: unanimous
        assert_eq!(result.sequence[1], b'C');
    }

    #[test]
    fn test_duplex_agreement() {
        let pos = ConsensusResult {
            sequence: b"ACGT".to_vec(),
            qualities: vec![40, 40, 40, 40],
            proportions: vec![1.0; 4],
            family_size: 3,
        };
        let neg = ConsensusResult {
            sequence: b"ACGT".to_vec(),
            qualities: vec![40, 40, 40, 40],
            proportions: vec![1.0; 4],
            family_size: 3,
        };
        let dcs = call_duplex_consensus_cpu(&pos, &neg);
        assert_eq!(dcs.sequence, b"ACGT");
    }

    #[test]
    fn test_duplex_disagreement() {
        let pos = ConsensusResult {
            sequence: b"ACGT".to_vec(),
            qualities: vec![40, 40, 40, 40],
            proportions: vec![1.0; 4],
            family_size: 3,
        };
        let neg = ConsensusResult {
            sequence: b"TCGT".to_vec(),
            qualities: vec![40, 40, 40, 40],
            proportions: vec![1.0; 4],
            family_size: 3,
        };
        let dcs = call_duplex_consensus_cpu(&pos, &neg);
        assert_eq!(dcs.sequence[0], b'N'); // A vs T -> N
        assert_eq!(dcs.sequence[1], b'C'); // agreement
    }
}
