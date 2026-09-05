#!/usr/bin/env python3
"""Decouple ANNOVAR database dir from install dir in annotate_variants.py."""
import ast, datetime, os, shutil, sys

TARGET = os.path.expanduser("~/pipelines/anneal/scripts/annotate_variants.py")
EDITS = [
    ("add DEFAULT_DB constant",
     'ANNOVAR_DB = os.path.join(ANNOVAR_DIR, "humandb")',
     'ANNOVAR_DB = os.path.join(ANNOVAR_DIR, "humandb")\n'
     '# Current databases live in references/humandb, NOT programs/annovar/humandb.\n'
     'DEFAULT_DB = os.path.join(\n'
     '    os.path.expanduser("~"), "references", "humandb")'),
    ("add --annovar-db argument",
     '''    ap.add_argument("--annovar-dir", default=ANNOVAR_DIR,
                    help="ANNOVAR installation directory")''',
     '''    ap.add_argument("--annovar-dir", default=ANNOVAR_DIR,
                    help="ANNOVAR installation directory (perl scripts)")
    ap.add_argument("--annovar-db", default=DEFAULT_DB,
                    help="ANNOVAR database directory (humandb with hg38_*.txt)")'''),
    ("run_annovar takes annovar_db explicitly",
     '''def run_annovar(vcf_in, out_prefix, annovar_dir):
    table_annovar = os.path.join(annovar_dir, "table_annovar.pl")
    annovar_db = os.path.join(annovar_dir, "humandb")
''',
     '''def run_annovar(vcf_in, out_prefix, annovar_dir, annovar_db):
    table_annovar = os.path.join(annovar_dir, "table_annovar.pl")
'''),
    ("main() passes annovar_db to run_annovar",
     "run_annovar(args.vcf, annovar_prefix, args.annovar_dir)",
     "run_annovar(args.vcf, annovar_prefix, args.annovar_dir, args.annovar_db)"),
]

def main():
    if not os.path.isfile(TARGET):
        sys.exit("Target not found: %s" % TARGET)
    src = open(TARGET).read()
    for desc, old, _ in EDITS:
        n = src.count(old)
        if n != 1:
            sys.exit("ABORT: edit '%s' matched %d times (expected 1)." % (desc, n))
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    bak = TARGET + ".bak_" + ts
    shutil.copy(TARGET, bak)
    print("Backup: %s" % bak)
    for desc, old, new in EDITS:
        src = src.replace(old, new, 1)
        print("Applied: %s" % desc)
    try:
        ast.parse(src)
    except SyntaxError as e:
        sys.exit("ABORT: syntax error after patch: %s" % e)
    open(TARGET, "w").write(src)
    print("Wrote patched %s\nSyntax OK." % TARGET)

if __name__ == "__main__":
    main()
