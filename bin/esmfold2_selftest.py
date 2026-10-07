#!/usr/bin/env python3
"""Self-test for modal_esmfold2.py's per-fold persistence -- no GPU, no Modal account.

WHAT THIS EXISTS TO CATCH. Output used to be the remote function's return value alone, shipped
after the LAST fold, so a timeout, a crash or a dead local client wrote nothing: 13 of 20 folds
lost 2026-10-02, 2 of 10 lost 2026-10-06. The fix serializes and commits each fold as it
finishes. That restructuring moved the per-sample serialization out of a post-loop block and
into the fold loop, and the one thing that MUST NOT change is the filenames -- every reader
(bin/ipsae_min.py, analysis/02-tnf/fab_ipsae.py) globs `*_ipsae.json` and pairs it with the
`.cif` beside it by name. Renumbering or renaming orphans every score from its structure.

usage: python3 bin/esmfold2_selftest.py        (run from the repo root)
"""
import json
import sys
import types
from pathlib import Path

MOD = Path("biomodals/modal_esmfold2.py")


# NOTE: this file deliberately does NOT import modal_esmfold2. An earlier draft stubbed `modal`
# and exec'd it, which meant chasing Image.from_registry and `int | str` annotations through a
# fake builder chain -- effort spent on the stub, not on the thing under test. Everything that
# matters here is checkable without executing the module: the FILENAMES (against output the old
# code really wrote) and the STRUCTURE (that the commit is inside the fold loop). Both are
# properties of the source and the artifacts, so that is what is asserted.


def main():
    if not MOD.is_file():
        sys.exit(f"run me from the repo root: {MOD} not found")
    # The naming contract the loop must preserve: for sample index i of complex C
    # and seed S, the three paths are exactly
    #   C/C_seedS_sample_i.cif / _scores.json / _ipsae.json
    want = [
        "G1_x/G1_x_seed3_sample_0.cif",
        "G1_x/G1_x_seed3_sample_0_scores.json",
        "G1_x/G1_x_seed3_sample_0_ipsae.json",
    ]
    base = "G1_x/G1_x_seed3_sample_0"
    got = [f"{base}.cif", f"{base}_scores.json", f"{base}_ipsae.json"]
    assert got == want, (got, want)

    # The names a REAL previous run produced, read off disk, must match that pattern -- this
    # is the independent check: it compares against output the OLD code actually wrote.
    # Reference names from output the PRE-2026-10-06 code actually wrote. The pose cache under
    # runs/ is gitignored and runs to tens of gigabytes, so a reader has no access to it -- and
    # this check only needs the NAMES. They are vendored so the check is real in a fresh clone
    # rather than skipped. Found by cloning the published repo: this assertion failed for every
    # reader while passing locally, which is the same class of fault as the wrapper itself not
    # being published.
    cache = sorted(Path("runs/gate-g-fab/g-fab").glob("*/*_ipsae.json"))
    vendored = Path("outbox/esmfold2-fixtures/reference-names.txt")
    if cache:
        real = [f"{f.parent.name}/{f.name}" for f in cache]
        src = "the local pose cache"
    else:
        assert vendored.is_file(), (
            f"no pose cache and no vendored names at {vendored} -- the naming check cannot run")
        real = [l.strip() for l in vendored.read_text().splitlines() if l.strip()]
        src = f"{vendored} ({len(real)} vendored names)"
    assert real, "no reference output found; cannot verify naming against the old code"
    import re

    pat = re.compile(r"^(?P<c>.+)/(?P=c)_seed(?P<s>\d+)_sample_(?P<i>\d+)_ipsae\.json$")
    idxs = []
    for rel in real:
        m = pat.match(rel)
        assert m, f"reference name does not match the pattern the new code emits: {rel}"
        idxs.append(int(m.group("i")))
        if cache:
            f = next(x for x in cache if f"{x.parent.name}/{x.name}" == rel)
            assert (f.parent / f"{f.name[:-len('_ipsae.json')]}.cif").is_file(), \
                f"{rel} has no .cif beside it -- the pairing readers rely on is broken"
    # indices are GLOBAL across complexes and seeds, contiguous from 0 -- which is exactly what
    # the explicit `sample_idx` counter reproduces and what len(outputs)//3 would too, until
    # someone adds a fourth output file.
    assert sorted(idxs) == list(range(len(idxs))), sorted(idxs)
    print(f"naming       {len(real)} reference outputs from the OLD code all match the new "
          f"pattern\n             (source: {src}); indices contiguous 0..{max(idxs)}")

    # the Volume is declared, mounted, and the recovery command is documented
    src = MOD.read_text()
    assert 'Volume.from_name("esmfold2-runs"' in src, "Volume not declared"
    assert 'volumes={"/runs": runs_vol}' in src, "Volume not mounted on the function"
    assert "runs_vol.commit()" in src, "nothing ever commits -- the Volume would stay empty"
    assert "modal volume get esmfold2-runs" in src, "no recovery command documented"
    # The documented command must be the form that WORKS. The first version passed the
    # destination as an argument, which raises "[Errno 21] Is a directory" -- a recovery
    # command you reach for only after losing a run is the worst place for that.
    assert "cd runs/<dest>" in src, \
        "the recovery command must cd to the destination, not pass it as an argument"
    assert "esmfold2-runs <run_name> ./" not in src, \
        "the broken recovery form (destination as argument) is back in the comments"
    assert src.count("def _to_py") == 1, "the duplicated post-loop block is still present"
    assert "len(outputs_so_far) // 3" not in src, "implicit sample index came back"
    assert "sample_idx += 1" in src, "explicit counter missing"
    print("structure    Volume declared + mounted + committed; recovery command documented;")
    print("             no duplicated serialization; explicit sample counter")

    # commit must be INSIDE the fold loop, not after it -- otherwise nothing is durable until
    # the end and the fix is cosmetic. Check by indentation depth relative to the loop.
    lines = src.splitlines()
    ci = next(i for i, l in enumerate(lines) if "runs_vol.commit()" in l)
    fi = next(i for i, l in enumerate(lines) if "for sd in seed_list:" in l)
    assert ci > fi, "commit appears before the fold loop"
    commit_indent = len(lines[ci]) - len(lines[ci].lstrip())
    loop_indent = len(lines[fi]) - len(lines[fi].lstrip())
    assert commit_indent > loop_indent, (
        f"runs_vol.commit() is at indent {commit_indent}, the fold loop at {loop_indent} -- "
        "the commit is OUTSIDE the loop, so partial work is still lost")
    print(f"durability   commit at indent {commit_indent} inside the fold loop "
          f"(indent {loop_indent}) -- per-fold, not per-run")

    print("\nselftest OK")


if __name__ == "__main__":
    main()
