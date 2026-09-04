#!/usr/bin/env python3
"""
patch_duplex_core.py

Applies targeted fixes to the anneal consensus engine so that DCS is
actually duplex.

Four defects addressed:

  1. families.rs   The CIGAR is accepted by build_family_tag() and never
                   written into the tag. Reads with different alignment
                   structure collapse into one family and are then
                   consensus-called by read index across misaligned
                   sequence. ConsensusCruncher's unique_tag() includes it.

  2. cpu.rs        call_duplex_consensus_cpu() gap-fills: when one strand
                   is N it emits the other strand's base. That is a
                   single-strand base labelled as duplex. ConsensusCruncher
                   emits N. Removed, and a per-strand Q30 floor added.

  3. correction.rs A rescued singleton inherits family_size = 1 + N from
                   call_duplex_consensus_cpu(), then that value is summed
                   again at DCS formation, so a 1-read strand can surface
                   as XW:i:41. Set to the true single-strand depth of 1.

  4. dcs.rs        No minimum per-strand support, and XW is a sum so a
                   20+1 pairing is indistinguishable from 20+20. Adds
                   config.min_reads_per_strand and emits XA/XB tags
                   carrying the two strands' depths separately.

Usage:
    python3 patch_duplex_core.py --repo ~/pipelines/anneal            # dry run
    python3 patch_duplex_core.py --repo ~/pipelines/anneal --apply    # write

Every edit is an exact-string replacement that must match exactly once.
If any anchor is missing or ambiguous the script aborts without writing.
Originals are copied to <file>.bak.<timestamp> before any change.
"""

import argparse
import os
import shutil
import sys
import time

# ---------------------------------------------------------------------------
# Edits: (relative path, description, old, new)
# ---------------------------------------------------------------------------

EDITS = []

# --- 1. families.rs: put the CIGAR back in the family key -------------------

EDITS.append((
    "src/grouping/families.rs",
    "family tag: add CIGAR field to format string",
    '"{}_{}_{}_{}_{}_{}_{}"',
    '"{}_{}_{}_{}_{}_{}_{}_{}"',
))

EDITS.append((
    "src/grouping/families.rs",
    "family tag: pass cigar into the format arguments",
    "barcode, chr1, start1, chr2, start2, strand, read_num",
    "barcode, chr1, start1, chr2, start2, cigar, strand, read_num",
))

EDITS.append((
    "src/grouping/families.rs",
    "test: complement tag must survive an embedded CIGAR",
    'let tag = "ACG+TTT_chr1_100_chr1_200_pos_R1";',
    'let tag = "ACG+TTT_chr1_100_chr1_200_100M_pos_R1";',
))

EDITS.append((
    "src/grouping/families.rs",
    "test: expected complement tag",
    'assert_eq!(complement, "TTT+ACG_chr1_100_chr1_200_neg_R2");',
    'assert_eq!(complement, "TTT+ACG_chr1_100_chr1_200_100M_neg_R2");',
))

# --- 2. cpu.rs: strict duplex, no gap-fill ---------------------------------

EDITS.append((
    "src/consensus/cpu.rs",
    "add duplex quality floor constant",
    "use crate::grouping::families::ReadData;",
    "use crate::grouping::families::ReadData;\n\n"
    "/// Minimum per-strand consensus quality for a duplex base to be called.\n"
    "/// Both strands must clear this independently.\n"
    "pub const DUPLEX_MIN_QUAL: u8 = 30;",
))

EDITS.append((
    "src/consensus/cpu.rs",
    "doc comment: remove the gap-fill contract",
    "///   - If one strand has N, use the other strand's base",
    "///   - If either strand is N, or either strand is below DUPLEX_MIN_QUAL,\n"
    "///     call N. There is no gap-fill: a base carried by only one strand\n"
    "///     is not a duplex observation.",
))

EDITS.append((
    "src/consensus/cpu.rs",
    "duplex consensus: delete gap-fill branches, require both strands",
    """        match (base_pos, base_neg) {
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
        }""",
    """        // Strict duplex. A base is emitted only when both strands
        // independently call the same non-N base at or above the quality
        // floor. Every other case -- disagreement, N on either strand,
        // low quality on either strand -- is N.
        let both_support = base_pos == base_neg
            && base_pos != b'N'
            && qual_pos >= DUPLEX_MIN_QUAL
            && qual_neg >= DUPLEX_MIN_QUAL;

        if both_support {
            let combined_qual = (qual_pos as u16 + qual_neg as u16).min(60) as u8;
            sequence.push(base_pos);
            qualities.push(combined_qual);
            proportions.push(1.0);
        } else {
            sequence.push(b'N');
            qualities.push(0);
            proportions.push(0.0);
        }""",
))

# --- 3. correction.rs: stop inflating rescued-singleton family size --------

EDITS.append((
    "src/singleton/correction.rs",
    "strategy 1: report true single-strand depth for rescued singleton",
    """                // Duplex-correct using the SSCS
                let corrected = call_duplex_consensus_cpu(
                    &singleton_consensus,
                    &complement_sscs.consensus,
                );""",
    """                // Duplex-correct using the SSCS
                let mut corrected = call_duplex_consensus_cpu(
                    &singleton_consensus,
                    &complement_sscs.consensus,
                );
                // call_duplex_consensus_cpu returns the SUM of both strands'
                // family sizes. This strand contributed exactly one read; if
                // the sum is kept it is summed a second time at DCS formation
                // and a 1-read strand surfaces as a large XW. Report the truth.
                corrected.family_size = 1;""",
))

EDITS.append((
    "src/singleton/correction.rs",
    "strategy 2: report true single-strand depth for singleton pair",
    """                    let corrected = call_duplex_consensus_cpu(&c1, &c2);""",
    """                    let mut corrected = call_duplex_consensus_cpu(&c1, &c2);
                    corrected.family_size = 1;""",
))

# --- 4. config.rs: per-strand minimum --------------------------------------

EDITS.append((
    "src/consensus/config.rs",
    "config: add min_reads_per_strand field",
    """    /// Enable singleton correction workflow.
    pub singleton_correction: bool,
}""",
    """    /// Enable singleton correction workflow.
    pub singleton_correction: bool,

    /// Minimum number of reads required on EACH strand before a DCS is
    /// formed. 1 reproduces the previous permissive behaviour; 2 or more
    /// enforces genuine two-sided duplex support.
    pub min_reads_per_strand: usize,
}""",
))

EDITS.append((
    "src/consensus/config.rs",
    "config: default min_reads_per_strand",
    """            singleton_correction: true,
        }""",
    """            singleton_correction: true,
            min_reads_per_strand: 1,
        }""",
))

# --- 5. dcs.rs: enforce per-strand support, carry both depths --------------

EDITS.append((
    "src/consensus/dcs.rs",
    "DcsRead: carry both strands' depths",
    """pub struct DcsRead {
    pub tag: String,
    pub consensus: crate::consensus::cpu::ConsensusResult,
    pub representative: ReadMetadata,
}""",
    """pub struct DcsRead {
    pub tag: String,
    pub consensus: crate::consensus::cpu::ConsensusResult,
    pub representative: ReadMetadata,
    /// Reads supporting the (+) strand SSCS of this duplex.
    pub fam_a: usize,
    /// Reads supporting the (-) strand SSCS of this duplex.
    pub fam_b: usize,
}""",
))

EDITS.append((
    "src/consensus/dcs.rs",
    "generate_dcs: use the config instead of ignoring it",
    "pub fn generate_dcs(sscs_reads: Vec<SscsRead>, _config: &ConsensusConfig) -> DcsResult {",
    "pub fn generate_dcs(sscs_reads: Vec<SscsRead>, config: &ConsensusConfig) -> DcsResult {",
))

EDITS.append((
    "src/consensus/dcs.rs",
    "generate_dcs: require min_reads_per_strand on both strands",
    """        // Mark both as paired
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
        });""",
    """        let a = by_tag.get(tag).unwrap();
        let b = by_tag.get(&comp).unwrap();
        let fam_a = a.consensus.family_size;
        let fam_b = b.consensus.family_size;

        // Strict duplex: both strands must independently carry enough reads.
        // A 20+1 pairing is one strand plus an anecdote, not duplex evidence.
        // Leave both tags unpaired so they fall through to sscs.singleton.
        if fam_a < config.min_reads_per_strand || fam_b < config.min_reads_per_strand {
            continue;
        }

        // Mark both as paired
        paired.insert(tag.clone());
        paired.insert(comp.clone());

        let duplex = call_duplex_consensus_cpu(&a.consensus, &b.consensus);
        let dcs_tag = format!("{}:{}+{}:{}", tag, fam_a, comp, fam_b);

        dcs_reads.push(DcsRead {
            tag: dcs_tag,
            consensus: duplex,
            representative: a.representative.clone(),
            fam_a,
            fam_b,
        });""",
))

# --- 6. pipeline.rs: emit XA/XB so per-strand depth is inspectable ---------

EDITS.append((
    "src/consensus/pipeline.rs",
    "write_dcs_sam: add XA/XB tag fields",
    '"{}\\t{}\\t{}\\t{}\\t{}\\t{}\\t{}\\t{}\\t{}\\t{}\\t{}\\tXV:Z:DCS\\tXW:i:{}",',
    '"{}\\t{}\\t{}\\t{}\\t{}\\t{}\\t{}\\t{}\\t{}\\t{}\\t{}\\tXV:Z:DCS\\tXW:i:{}\\tXA:i:{}\\tXB:i:{}",',
))

EDITS.append((
    "src/consensus/pipeline.rs",
    "write_dcs_sam: supply XA/XB values",
    """        dcs.consensus.family_size,
    )?;""",
    """        dcs.consensus.family_size,
        dcs.fam_a,
        dcs.fam_b,
    )?;""",
))

# --- 7. main.rs: wire min_reads_per_strand --------------------------------

EDITS.append((
    "src/main.rs",
    "run_stage2: set min_reads_per_strand from ANNEAL_MIN_READS_PER_STRAND",
    """        use_gpu: gpu_available,
        gpu_device,
        singleton_correction,
    };""",
    """        use_gpu: gpu_available,
        gpu_device,
        singleton_correction,
        min_reads_per_strand: std::env::var("ANNEAL_MIN_READS_PER_STRAND")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .filter(|v| *v >= 1)
            .unwrap_or(1),
    };""",
))


# ---------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True,
                    help="Path to the anneal working tree")
    ap.add_argument("--apply", action="store_true",
                    help="Write changes. Without this the script only reports.")
    args = ap.parse_args()

    repo = os.path.abspath(os.path.expanduser(args.repo))
    if not os.path.isdir(os.path.join(repo, "src")):
        sys.exit("ERROR: %s does not look like the anneal tree (no src/)" % repo)

    # Load every file once.
    contents = {}
    for rel, _desc, _old, _new in EDITS:
        if rel in contents:
            continue
        path = os.path.join(repo, rel)
        if not os.path.isfile(path):
            sys.exit("ERROR: missing file %s" % path)
        with open(path, "r", encoding="utf-8") as fh:
            contents[rel] = fh.read()

    # Verify and stage every edit before writing anything.
    failures = []
    for rel, desc, old, new in EDITS:
        text = contents[rel]
        n = text.count(old)
        if n == 1:
            contents[rel] = text.replace(old, new, 1)
            print("  OK        %-28s %s" % (rel, desc))
        elif n == 0:
            if text.count(new) >= 1:
                print("  SKIP      %-28s %s (already applied)" % (rel, desc))
            else:
                print("  NOT FOUND %-28s %s" % (rel, desc))
                failures.append((rel, desc, "anchor not present"))
        else:
            print("  AMBIGUOUS %-28s %s (%d matches)" % (rel, desc, n))
            failures.append((rel, desc, "%d matches" % n))

    print("")
    if failures:
        print("ABORTED. %d edit(s) could not be applied cleanly:" % len(failures))
        for rel, desc, why in failures:
            print("    %s :: %s :: %s" % (rel, desc, why))
        print("")
        print("The working tree has NOT been modified. This usually means the")
        print("local source has diverged from the version this patch targets.")
        print("Send me the affected file and I will re-anchor the edit.")
        sys.exit(1)

    if not args.apply:
        print("Dry run only. All %d edits verified. Re-run with --apply to write."
              % len(EDITS))
        return

    stamp = time.strftime("%Y%m%d-%H%M%S")
    for rel, text in contents.items():
        path = os.path.join(repo, rel)
        backup = "%s.bak.%s" % (path, stamp)
        shutil.copy2(path, backup)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
        print("  wrote %s   (backup: %s)" % (rel, os.path.basename(backup)))

    print("")
    print("Done. %d files changed." % len(contents))
    print("")
    print("Next:")
    print("    cd %s" % repo)
    print("    CARGO_TARGET_DIR=target_cpu cargo build --release")
    print("    CARGO_TARGET_DIR=target_cpu cargo test --release")


if __name__ == "__main__":
    main()
