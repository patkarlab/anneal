#!/usr/bin/env python3
"""
patch_vep_env.py

Restores the VEP environment sanitization to scripts/annotate_variants.py.

History: the June 2026 annotation work solved this in bin/annotate_anneal.py,
which ran from a standalone PBS array with no conda env active. The fix never
reached annotate_variants.py, which has no way to pass an environment to its
subprocesses. Stage 3 now runs inside run_pipeline.sh, which calls
activate_conda -> `conda activate anneal`, putting envs/anneal/bin on PATH. Its
perl then shadows VEP's own and every VEP module load fails with:

    Compilation failed in require at .../envs/anneal/lib/perl5/core_perl/base.pm

Two changes:

  1. run() takes an optional env, passed through to subprocess.run.

  2. run_vep() reverts to `conda run -n vep vep` (calling the binary directly
     hits the same perl shadowing, plus `conda run -n` alone resolves to
     ~/.conda/envs and fails) and builds a sanitized environment: any conda env
     bin that would shadow VEP's perl is dropped from PATH, and the leaked Perl
     library variables are cleared.

The strip list covers targeted-seq (the original culprit) and anneal (the
current one). Add to VEP_PATH_STRIP if another env ever shadows it.
"""

import shutil
import sys
import time

P = "scripts/annotate_variants.py"

RUN_OLD = '''def run(cmd, desc=None, shell=False):
    if desc:
        log.info("%s", desc)
    cmd_str = cmd if shell else " ".join(cmd)
    log.info("  cmd: %s", cmd_str)
    result = subprocess.run(cmd, capture_output=True, text=True, shell=shell)'''

RUN_NEW = '''# Conda env bin directories whose perl would shadow VEP's own. VEP loads its
# modules relative to whichever perl is found first on PATH, so any of these
# appearing ahead of the vep env breaks every module load.
VEP_PATH_STRIP = ("envs/targeted-seq/bin", "envs/anneal/bin")

# Perl library variables that leak in from an activated env and misdirect VEP.
VEP_PERL_VARS = ("PERL5LIB", "PERL_LOCAL_LIB_ROOT", "PERL_MM_OPT", "PERL_MB_OPT")


def vep_environment():
    """A copy of the current environment safe for VEP to run in."""
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join(
        p for p in env.get("PATH", "").split(os.pathsep)
        if not any(s in p for s in VEP_PATH_STRIP))
    for v in VEP_PERL_VARS:
        env.pop(v, None)
    return env


def run(cmd, desc=None, shell=False, env=None):
    if desc:
        log.info("%s", desc)
    cmd_str = cmd if shell else " ".join(cmd)
    log.info("  cmd: %s", cmd_str)
    result = subprocess.run(cmd, capture_output=True, text=True, shell=shell,
                            env=env)'''

VEP_OLD = '''    cmd = [
        # Called directly rather than through `conda run`: a subprocess does
        # not inherit conda's env-path config, so `conda run -n vep` resolves to
        # ~/.conda/envs/vep and fails with EnvironmentLocationNotFound. Override
        # with the VEP_BIN environment variable.
        os.environ.get("VEP_BIN",
                       os.path.expanduser("~/miniconda3/envs/vep/bin/vep")),'''

VEP_NEW = '''    cmd = [
        # Must go through `conda run` so VEP gets its own perl. Calling the
        # binary directly leaves the active env's perl on PATH, which shadows
        # VEP's modules. `-p <prefix>` rather than `-n <name>` because a
        # subprocess does not inherit conda's env-path config and `-n vep`
        # resolves to ~/.conda/envs/vep, which does not exist.
        "conda", "run", "-p",
        os.environ.get("VEP_PREFIX",
                       os.path.expanduser("~/miniconda3/envs/vep")),
        "vep",'''


def main():
    try:
        s = open(P).read()
    except FileNotFoundError:
        sys.exit(f"FATAL: run from the anneal root; {P} not found")

    if "def vep_environment" in s:
        sys.exit("Already patched. Nothing to do.")

    for name, old in (("run()", RUN_OLD), ("run_vep() command", VEP_OLD)):
        n = s.count(old)
        if n != 1:
            sys.stderr.write(f"FATAL: {name} anchor found {n} times, expected 1\n")
            sys.stderr.write(old[:200] + "\n")
            sys.exit(1)
        print(f"anchor ok: {name}")

    s = s.replace(RUN_OLD, RUN_NEW, 1)
    s = s.replace(VEP_OLD, VEP_NEW, 1)

    # route the VEP call through the sanitized environment
    old_call = 'return run(cmd, desc="Running VEP on " + os.path.basename(vcf_in))'
    if old_call in s:
        s = s.replace(old_call,
                      'return run(cmd, desc="Running VEP on " + os.path.basename(vcf_in),\n'
                      '               env=vep_environment())', 1)
        print("anchor ok: VEP run() call")
    else:
        print("NOTE: could not find the VEP run() call; add env=vep_environment() by hand")

    backup = f"{P}.bak_{time.strftime('%Y%m%d_%H%M%S')}"
    shutil.copy2(P, backup)
    open(P, "w").write(s)
    print(f"\npatched {P}\nbackup  {backup}")


if __name__ == "__main__":
    main()
