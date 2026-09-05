#!/usr/bin/env python3
"""
consensus_callers.py — Blind three-caller consensus for the anneal MRD pipeline.

Ingests VarDict, LoFreq (baq_off), and Pisces (gVCF) variant calls on the SAME
panel BAM, normalizes them to a common schema, and emits a consensus table with a
caller-count column (1/2/3 = how many callers agree on each variant). Optionally
runs ANNOVAR and merges the gene annotation back onto the consensus.

BLIND BY DESIGN: this script never takes a diagnosis-variant list. It reconciles
callers and counts agreement with no knowledge of which variants are "expected".
Recovery scoring against a truth set is a SEPARATE downstream step.

Representation handling:
  - SNVs match on exact (chrom, pos, ref, alt).
  - Indels/complex match on a locus window (chrom, pos within +/- WINDOW), so the
    three callers' different indel anchorings (e.g. NPM1) count as concordant.
    Window-matched indels are flagged so the merge is auditable.

Strand data differs by caller (LoFreq DP4, VarDict VARBIAS/ALD, Pisces AD+SB), so
per-caller strand info is recorded where available but NOT forced into a uniform
filter. Caller-count is the primary consensus axis.
"""

import argparse
import os
import re
import sys
from collections import defaultdict


# ----------------------------- VCF parsing -----------------------------

def _info_dict(info):
    d = {}
    for kv in info.split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k] = v
        elif kv:
            d[kv] = True
    return d


def _fmt_dict(fmt, sample):
    keys = fmt.split(':')
    vals = sample.split(':')
    return dict(zip(keys, vals))


def _is_indel(ref, alt):
    return len(ref) != len(alt)


def parse_vardict(path):
    """VarDict VCF -> list of variant dicts. AF from FORMAT or INFO; VARBIAS=fwd:rev alt strands."""
    out = []
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 8:
                continue
            chrom, pos, _id, ref, alt, qual, filt, info = f[:8]
            if alt == '.' or ',' in alt:  # skip non-calls and multiallelic (rare; split upstream if needed)
                if ',' in alt:
                    pass  # could split; for now take first allele path below
            info_d = _info_dict(info)
            af = None
            alt_fwd = alt_rev = None
            # FORMAT/SAMPLE if present
            if len(f) >= 10:
                fmt_d = _fmt_dict(f[8], f[9])
                if 'AF' in fmt_d:
                    try: af = float(fmt_d['AF'])
                    except ValueError: pass
                # ALD = alt fwd,rev (VarDict)
                if 'ALD' in fmt_d and ',' in fmt_d['ALD']:
                    a, b = fmt_d['ALD'].split(',')[:2]
                    try: alt_fwd, alt_rev = int(a), int(b)
                    except ValueError: pass
            if af is None and 'AF' in info_d:
                try: af = float(info_d['AF'])
                except (ValueError, TypeError): pass
            # VARBIAS in INFO = "reffwd:refrev:altfwd:altrev" or "fwd:rev" depending on version
            if (alt_fwd is None) and 'VARBIAS' in info_d:
                parts = re.split('[:,]', str(info_d['VARBIAS']))
                if len(parts) >= 2:
                    try:
                        alt_fwd, alt_rev = int(parts[-2]), int(parts[-1])
                    except ValueError:
                        pass
            out.append(dict(chrom=chrom, pos=int(pos), ref=ref, alt=alt,
                            af=af, alt_fwd=alt_fwd, alt_rev=alt_rev,
                            filter=filt, caller='vardict'))
    return out


def parse_lofreq(path):
    """LoFreq VCF -> variant dicts. INFO AF; DP4=ref-fwd,ref-rev,alt-fwd,alt-rev."""
    out = []
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 8:
                continue
            chrom, pos, _id, ref, alt, qual, filt, info = f[:8]
            if alt == '.':
                continue
            info_d = _info_dict(info)
            af = None
            if 'AF' in info_d:
                try: af = float(info_d['AF'])
                except (ValueError, TypeError): pass
            alt_fwd = alt_rev = None
            if 'DP4' in info_d:
                p = info_d['DP4'].split(',')
                if len(p) == 4:
                    try: alt_fwd, alt_rev = int(p[2]), int(p[3])
                    except ValueError: pass
            out.append(dict(chrom=chrom, pos=int(pos), ref=ref, alt=alt,
                            af=af, alt_fwd=alt_fwd, alt_rev=alt_rev,
                            filter=filt, caller='lofreq'))
    return out


def parse_pisces(path):
    """Pisces gVCF -> variant dicts (ALT != '.'). FORMAT: AD=ref,alt ; VF=freq ; SB=bias score.
    Pisces gives no fwd/rev split, so alt_fwd/alt_rev stay None; we keep alt depth + SB."""
    out = []
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 8:
                continue
            chrom, pos, _id, ref, alt, qual, filt, info = f[:8]
            if alt == '.' or alt == '<M>':  # gVCF reference / no-variant block
                continue
            af = None; alt_depth = None; sb = None
            if len(f) >= 10:
                fmt_d = _fmt_dict(f[8], f[9])
                if 'VF' in fmt_d:
                    try: af = float(fmt_d['VF'])
                    except ValueError: pass
                if 'AD' in fmt_d and ',' in fmt_d['AD']:
                    try: alt_depth = int(fmt_d['AD'].split(',')[1])
                    except (ValueError, IndexError): pass
                if 'SB' in fmt_d:
                    try: sb = float(fmt_d['SB'])
                    except ValueError: pass
            out.append(dict(chrom=chrom, pos=int(pos), ref=ref, alt=alt,
                            af=af, alt_fwd=None, alt_rev=None, alt_depth=alt_depth,
                            sb=sb, filter=filt, caller='pisces'))
    return out


# ----------------------------- Matching -----------------------------

def variant_key(v):
    """Exact key for SNVs/MNVs."""
    return (v['chrom'], v['pos'], v['ref'].upper(), v['alt'].upper())


def build_consensus(by_caller, indel_window=5, pass_only_pisces=False):
    """
    by_caller: dict caller_name -> list of variant dicts.
    Returns list of consensus records, each with per-caller presence + caller_count.

    SNVs: exact (chrom,pos,ref,alt) match.
    Indels: a variant from caller B is merged into an existing indel record from
            caller A if same chrom and |pos_A - pos_B| <= indel_window and both indels.
            (Different anchoring of the same event -> concordant.)
    """
    callers = list(by_caller.keys())

    # 1) Index SNVs by exact key; collect indels separately for window matching.
    snv_records = {}   # key -> record
    indel_records = [] # list of records (each may accumulate multiple callers)

    def new_record(v, is_indel):
        rec = dict(chrom=v['chrom'], pos=v['pos'], ref=v['ref'], alt=v['alt'],
                   is_indel=is_indel, window_matched=False,
                   callers=set(), per_caller={})
        return rec

    def add_to_record(rec, v):
        c = v['caller']
        rec['callers'].add(c)
        rec['per_caller'][c] = v

    for caller in callers:
        for v in by_caller[caller]:
            if pass_only_pisces and caller == 'pisces' and v.get('filter') != 'PASS':
                continue
            indel = _is_indel(v['ref'], v['alt'])
            if not indel:
                k = variant_key(v)
                if k not in snv_records:
                    snv_records[k] = new_record(v, False)
                add_to_record(snv_records[k], v)
            else:
                # try to merge into an existing indel record within window
                merged = False
                for rec in indel_records:
                    if rec['chrom'] != v['chrom']:
                        continue
                    if abs(rec['pos'] - v['pos']) <= indel_window:
                        add_to_record(rec, v)
                        if v['pos'] != rec['pos'] or v['alt'] != rec['alt']:
                            rec['window_matched'] = True
                        merged = True
                        break
                if not merged:
                    indel_records.append(new_record(v, True))
                    add_to_record(indel_records[-1], v)

    records = list(snv_records.values()) + indel_records
    for rec in records:
        rec['caller_count'] = len(rec['callers'])
    records.sort(key=lambda r: (r['chrom'], r['pos']))
    return records, callers


# ----------------------------- Output -----------------------------

def af_of(rec, caller):
    v = rec['per_caller'].get(caller)
    if v and v.get('af') is not None:
        return f"{v['af']:.6f}"
    return '.'

def strand_of(rec, caller):
    v = rec['per_caller'].get(caller)
    if not v:
        return '.'
    if v.get('alt_fwd') is not None and v.get('alt_rev') is not None:
        return f"{v['alt_fwd']}:{v['alt_rev']}"
    if v.get('alt_depth') is not None:  # Pisces: alt depth + SB
        sb = v.get('sb')
        return f"AD={v['alt_depth']}" + (f";SB={sb}" if sb is not None else "")
    return '.'

def write_consensus_tsv(records, callers, path):
    cols = ['chrom', 'pos', 'ref', 'alt', 'is_indel', 'window_matched', 'caller_count']
    for c in callers:
        cols += [f'{c}_AF', f'{c}_strand']
    cols += ['callers']
    with open(path, 'w') as out:
        out.write('\t'.join(cols) + '\n')
        for r in records:
            row = [r['chrom'], str(r['pos']), r['ref'], r['alt'],
                   '1' if r['is_indel'] else '0',
                   '1' if r['window_matched'] else '0',
                   str(r['caller_count'])]
            for c in callers:
                row += [af_of(r, c), strand_of(r, c)]
            row += [','.join(sorted(r['callers']))]
            out.write('\t'.join(row) + '\n')


def write_avinput(records, path):
    """ANNOVAR avinput: chrom start end ref alt  (1-based; convention used by table_annovar).
    For indels we emit ANNOVAR-style: insertion ref='-', deletion alt='-'. We pass the
    VCF-style anchored alleles through annovar's -vcfinput instead, so here we emit a
    minimal VCF for annotation to avoid hand-converting indel representations."""
    with open(path, 'w') as out:
        out.write('##fileformat=VCFv4.1\n')
        out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for r in records:
            out.write('\t'.join([r['chrom'], str(r['pos']), '.', r['ref'], r['alt'],
                                 '.', '.', f"CC={r['caller_count']}"]) + '\n')


def merge_annovar(records, multianno_path):
    """Merge ANNOVAR multianno gene annotation back onto consensus records by (chrom,pos,ref,alt)."""
    if not multianno_path or not os.path.exists(multianno_path):
        return None
    ann = {}
    with open(multianno_path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        # Identify columns
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.rstrip('\n').split('\t')
            try:
                key = (f[idx['Chr']], f[idx['Start']], f[idx['Ref']], f[idx['Alt']])
            except (KeyError, IndexError):
                continue
            ann[key] = {
                'Func': f[idx.get('Func.refGene', -1)] if 'Func.refGene' in idx else '.',
                'Gene': f[idx.get('Gene.refGene', -1)] if 'Gene.refGene' in idx else '.',
                'ExonicFunc': f[idx.get('ExonicFunc.refGene', -1)] if 'ExonicFunc.refGene' in idx else '.',
                'AAChange': f[idx.get('AAChange.refGene', -1)] if 'AAChange.refGene' in idx else '.',
            }
    return ann


def write_annotated(records, callers, ann, path):
    cols = ['chrom', 'pos', 'ref', 'alt', 'caller_count',
            'Func.refGene', 'Gene.refGene', 'ExonicFunc.refGene', 'AAChange.refGene']
    for c in callers:
        cols += [f'{c}_AF']
    cols += ['is_indel', 'window_matched', 'callers']
    with open(path, 'w') as out:
        out.write('\t'.join(cols) + '\n')
        for r in records:
            # ANNOVAR may have normalized indel coords; try exact, else annotate as '.'
            a = None
            if ann is not None:
                a = ann.get((r['chrom'], str(r['pos']), r['ref'], r['alt']))
            a = a or {'Func': '.', 'Gene': '.', 'ExonicFunc': '.', 'AAChange': '.'}
            row = [r['chrom'], str(r['pos']), r['ref'], r['alt'], str(r['caller_count']),
                   a['Func'], a['Gene'], a['ExonicFunc'], a['AAChange']]
            for c in callers:
                row += [af_of(r, c)]
            row += ['1' if r['is_indel'] else '0',
                    '1' if r['window_matched'] else '0',
                    ','.join(sorted(r['callers']))]
            out.write('\t'.join(row) + '\n')


# ----------------------------- CLI -----------------------------

def main():
    ap = argparse.ArgumentParser(description="Blind three-caller consensus with caller-count + ANNOVAR.")
    ap.add_argument('--vardict', help='VarDict VCF')
    ap.add_argument('--lofreq', help='LoFreq baq_off VCF')
    ap.add_argument('--pisces', help='Pisces gVCF')
    ap.add_argument('--out-prefix', required=True, help='Output prefix')
    ap.add_argument('--indel-window', type=int, default=5,
                    help='bp window for indel anchoring reconciliation (default 5)')
    ap.add_argument('--pisces-pass-only', action='store_true',
                    help='Only count Pisces variants with FILTER=PASS')
    ap.add_argument('--multianno', help='ANNOVAR hg38_multianno.txt to merge (optional)')
    args = ap.parse_args()

    by_caller = {}
    if args.vardict: by_caller['vardict'] = parse_vardict(args.vardict)
    if args.lofreq:  by_caller['lofreq']  = parse_lofreq(args.lofreq)
    if args.pisces:  by_caller['pisces']  = parse_pisces(args.pisces)
    if not by_caller:
        sys.exit("ERROR: provide at least one caller VCF")

    for c, vs in by_caller.items():
        sys.stderr.write(f"{c}: {len(vs)} variant records\n")

    records, callers = build_consensus(by_caller,
                                       indel_window=args.indel_window,
                                       pass_only_pisces=args.pisces_pass_only)

    # consensus TSV
    tsv = f"{args.out_prefix}.consensus.tsv"
    write_consensus_tsv(records, callers, tsv)
    # avinput-style VCF for ANNOVAR
    annvcf = f"{args.out_prefix}.for_annovar.vcf"
    write_avinput(records, annvcf)

    # counts by caller_count
    counts = defaultdict(int)
    for r in records:
        counts[r['caller_count']] += 1
    sys.stderr.write("\n=== consensus summary ===\n")
    sys.stderr.write(f"total unique variants: {len(records)}\n")
    for k in sorted(counts, reverse=True):
        sys.stderr.write(f"  called by {k} caller(s): {counts[k]}\n")
    sys.stderr.write(f"\nconsensus table: {tsv}\n")
    sys.stderr.write(f"annovar input:   {annvcf}\n")

    # optional annotation merge
    if args.multianno:
        ann = merge_annovar(records, args.multianno)
        annotated = f"{args.out_prefix}.consensus.annotated.tsv"
        write_annotated(records, callers, ann, annotated)
        sys.stderr.write(f"annotated table: {annotated}\n")


if __name__ == '__main__':
    main()
