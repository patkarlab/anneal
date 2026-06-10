//! Read family construction: group aligned reads by UMI + genomic position.
//!
//! A "read family" shares: UMI barcode, read chromosome, read start position,
//! mate chromosome, mate start position, CIGAR string, strand orientation,
//! and read number. This composite key ensures reads from the same PCR
//! duplicate cluster are grouped together.

use ahash::AHashMap;
use log::debug;

/// Compact representation of a read for consensus calling.
/// Avoids holding full pysam-equivalent objects in memory.
#[derive(Debug, Clone)]
pub struct ReadData {
    /// Original query name (for tracking)
    pub qname: String,
    /// Consensus family tag
    pub tag: String,
    /// DNA sequence as bytes (ACGTN)
    pub sequence: Vec<u8>,
    /// Phred quality scores (raw, not +33 encoded)
    pub qualities: Vec<u8>,
    /// SAM flag
    pub flag: u16,
    /// Reference sequence name (chromosome)
    pub rname: String,
    /// Reference sequence ID (chromosome index)
    pub ref_id: i32,
    /// 0-based leftmost mapping position
    pub pos: i64,
    /// Mapping quality
    pub mapq: u8,
    /// CIGAR string (as text for tag construction)
    pub cigar_str: String,
    /// Mate reference name
    pub mate_rname: String,
    /// Mate reference ID
    pub mate_ref_id: i32,
    /// Mate position
    pub mate_pos: i64,
    /// Template length
    pub template_len: i64,
}

/// A family of reads sharing the same UMI + genomic coordinates.
#[derive(Debug)]
pub struct ReadFamily {
    /// The consensus tag identifying this family
    pub tag: String,
    /// All reads in this family
    pub reads: Vec<ReadData>,
}

impl ReadFamily {
    pub fn family_size(&self) -> usize {
        self.reads.len()
    }

    pub fn is_singleton(&self) -> bool {
        self.reads.len() == 1
    }
}

/// Determines strand orientation from SAM flags.
/// Returns "pos" or "neg" based on read pair orientation.
pub fn determine_strand(flag: u16) -> &'static str {
    // Properly paired reads: 99/147 = pos strand, 83/163 = neg strand
    match flag {
        99 | 147 => "pos",
        83 | 163 => "neg",
        // For edge cases (translocations, etc.), use bit flags
        _ => {
            if flag & 0x10 != 0 {
                "neg"
            } else {
                "pos"
            }
        }
    }
}

/// Determines read number from SAM flags.
pub fn determine_read_number(flag: u16) -> &'static str {
    if flag & 0x40 != 0 {
        "R1"
    } else if flag & 0x80 != 0 {
        "R2"
    } else {
        "R0"
    }
}

/// Constructs the consensus family tag from read attributes.
///
/// Format: [Barcode]_[ReadChr]_[ReadStart]_[MateChr]_[MateStart]_[CigarOrdered]_[Strand]_[ReadNum]
///
/// The barcode is extracted from the query name (appended during barcode extraction).
/// CIGAR strings are ordered by strand and read number for consistent pairing.
pub fn build_family_tag(
    qname: &str,
    ref_name: &str,
    pos: i64,
    mate_ref_name: &str,
    mate_pos: i64,
    cigar: &str,
    flag: u16,
    barcode_delim: char,
) -> String {
    // Extract barcode from query name: "READNAME|UMI1+UMI2" -> "UMI1+UMI2"
    let barcode = qname
        .rsplit(barcode_delim)
        .next()
        .unwrap_or("UNKNOWN");

    let strand = determine_strand(flag);
    let read_num = determine_read_number(flag);

    // Order coordinates: always put lower chr/pos first for consistent grouping
    let (chr1, start1, chr2, start2) = if ref_name < mate_ref_name
        || (ref_name == mate_ref_name && pos <= mate_pos)
    {
        (ref_name, pos, mate_ref_name, mate_pos)
    } else {
        (mate_ref_name, mate_pos, ref_name, pos)
    };

    format!(
        "{}_{}_{}_{}_{}_{}_{}", 
        barcode, chr1, start1, chr2, start2, strand, read_num
    )
}

/// Builds the duplex complement tag by swapping the barcode halves
/// and flipping the strand orientation.
///
/// For a tag "ACG+TTT_chr1_100_chr1_200_98M_pos",
/// the duplex complement is "TTT+ACG_chr1_100_chr1_200_98M_neg".
pub fn build_duplex_complement_tag(tag: &str) -> Option<String> {
    // Tag format: BARCODE_chr1_pos1_chr2_pos2_cigar_strand
    let parts: Vec<&str> = tag.splitn(2, '_').collect();
    if parts.len() < 2 {
        return None;
    }

    let barcode = parts[0];
    let rest = parts[1];

    // Swap barcode halves: "UMI1+UMI2" -> "UMI2+UMI1"
    let barcode_parts: Vec<&str> = barcode.split('+').collect();
    if barcode_parts.len() != 2 {
        return None;
    }
    let swapped_barcode = format!("{}+{}", barcode_parts[1], barcode_parts[0]);

    // Duplex pairing geometry:
    //   pos_R1 <-> neg_R2 (both cover left end of insert)
    //   pos_R2 <-> neg_R1 (both cover right end of insert)
    // So we flip strand AND swap read number.
    let mut flipped = rest.to_string();

    // Flip strand
    if flipped.contains("_pos_") {
        flipped = flipped.replacen("_pos_", "_neg_", 1);
    } else if flipped.contains("_neg_") {
        flipped = flipped.replacen("_neg_", "_pos_", 1);
    } else {
        return None;
    }

    // Swap read number
    if flipped.ends_with("_R1") {
        flipped = format!("{}_R2", &flipped[..flipped.len() - 3]);
    } else if flipped.ends_with("_R2") {
        flipped = format!("{}_R1", &flipped[..flipped.len() - 3]);
    }

    Some(format!("{}_{}", swapped_barcode, flipped))
}

/// Group reads from a genomic region into families.
///
/// Takes an iterator of ReadData and returns families keyed by their tag.
/// This is the core grouping operation that determines consensus families.
pub fn group_reads_into_families(
    reads: impl Iterator<Item = ReadData>,
    max_family_size: usize,
) -> AHashMap<String, ReadFamily> {
    let mut families: AHashMap<String, Vec<ReadData>> = AHashMap::new();

    for read in reads {
        let tag = read.tag.clone();
        families.entry(tag).or_default().push(read);
    }

    // Convert to ReadFamily, applying max family size cap
    families
        .into_iter()
        .map(|(tag, mut reads)| {
            if reads.len() > max_family_size {
                debug!(
                    "Family {} has {} reads, downsampling to {}",
                    tag,
                    reads.len(),
                    max_family_size
                );
                reads.truncate(max_family_size);
            }
            (
                tag.clone(),
                ReadFamily { tag, reads },
            )
        })
        .collect()
}

/// Separate families into SSCS-eligible (size >= min) and singletons.
pub fn partition_families(
    families: AHashMap<String, ReadFamily>,
    min_family_size: usize,
) -> (Vec<ReadFamily>, Vec<ReadFamily>) {
    let mut sscs_families = Vec::new();
    let mut singletons = Vec::new();

    for (_, family) in families {
        if family.family_size() >= min_family_size && family.family_size() > 1 {
            sscs_families.push(family);
        } else {
            singletons.push(family);
        }
    }

    (sscs_families, singletons)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_duplex_complement_tag() {
        let tag = "ACG+TTT_chr1_100_chr1_200_pos_R1";
        let complement = build_duplex_complement_tag(tag).unwrap();
        assert_eq!(complement, "TTT+ACG_chr1_100_chr1_200_neg_R2");
    }

    #[test]
    fn test_strand_determination() {
        assert_eq!(determine_strand(99), "pos");
        assert_eq!(determine_strand(147), "pos");
        assert_eq!(determine_strand(83), "neg");
        assert_eq!(determine_strand(163), "neg");
    }
}
