#!/usr/bin/env python3
"""Measure NPM1 and CEBPA insertion VAF across the G/H/I/J dilution rungs.
Indel-aware (pysam pileup), so it sees what count_coverage cannot.
Tells us whether these known drivers DILUTE (real, usable controls) or are flat."""
import pysam, glob

def count_insertion(bam_path, chrom, pos1, inserted_seq, min_bq=30):
    p0 = pos1 - 1
    bam = pysam.AlignmentFile(bam_path)
    alt = depth = 0
    for col in bam.pileup(chrom, p0, p0+1, truncate=True, max_depth=500000,
                          min_base_quality=0, stepper="nofilter"):
        for r in col.pileups:
            if r.query_position is None:
                # still counts toward depth if it's a covering read
                depth += 1; continue
            # base-quality gate on the anchor base
            bq = r.alignment.query_qualities[r.query_position]
            if bq < min_bq:
                continue
            depth += 1
            if r.indel == len(inserted_seq):
                q = r.alignment.query_sequence
                ins = q[r.query_position+1 : r.query_position+1+r.indel]
                if ins.upper() == inserted_seq.upper():
                    alt += 1
    return alt, depth

drivers = [
    ("NPM1 type A  chr5:171410539 C>CTCTG", "chr5",  171410539, "TCTG"),
    ("CEBPA        chr19:33302346 C>CG",    "chr19", 33302346,  "G"),
]
print(f"{'driver':40} {'rung':4} {'altIns':>7} {'depth':>7} {'VAF%':>8}  fold-vs-prev")
for label, chrom, pos1, ins in drivers:
    prev=None
    for R in "GHIJ":
        b=glob.glob(f"results_dilution_gpu/DIL-A-{R}-Duplex/consensus/DIL-A-{R}-Duplex.sscs.sc.sorted.bam")[0]
        a,d = count_insertion(b, chrom, pos1, ins)
        vaf = 100*a/d if d else 0
        fold = "" if prev is None else (f"{prev/vaf:.1f}x" if vaf>0 else "->floor")
        print(f"{label:40} {R:4} {a:7d} {d:7d} {vaf:8.4f}  {fold}")
        prev = vaf
    print()
