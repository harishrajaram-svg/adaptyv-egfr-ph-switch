#!/usr/bin/env python3
"""Regression tests for the two fail-OPEN faults PK found on 2026-10-05.

    python3 test_failclosed.py        # exits non-zero if either fault returns

Both faults shared a shape: the scorer returned a NUMBER where it had no valid
measurement. Neither changed a submitted value (verified: all 6840 cached outputs
are single-pair/two-direction, so neither path was ever taken), but "latent" is
not "absent", and a scorer that can invent a number cannot be the reproducible
path the submission rests on.

  FAULT 1 -- stale output survives a failed run.
    bin/ipsae_min.score() ran the reference, then read <stem>_10_10.txt. It
    checked neither the return code nor whether that file was written by THIS
    run. A crashed re-score therefore returned the PREVIOUS run's number.

  FAULT 2 -- one direction accepted as a minimum.
    ipSAE_min is the minimum over an interface's two reciprocal alignment
    directions. A file holding only one direction was accepted, and that single
    value reported as if it were the min. master_rank.py did this explicitly
    (`return vals[0] if vals else None`) and it writes the canonical file.
"""
import importlib.util, json, shutil, sys, tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent

def load(name):
    path = REPO / "bin" / f"{name}.py"
    if not path.exists(): path = HERE / "vendor" / "bin" / f"{name}.py"
    if not path.exists(): sys.exit(f"cannot find {name}.py")
    spec = importlib.util.spec_from_file_location(f"_t_{name}", path)
    m = importlib.util.module_from_spec(spec); sys.modules[spec.name] = m
    spec.loader.exec_module(m); return m

LIVE, SEEDAGG, CANON = load("ipsae_min"), load("instrument_v2"), load("master_rank")
HEADER = ("Chn1 Chn2 PAE Dist Type ipSAE ipSAE_d0chn ipSAE_d0dom ipTM_af ipTM_d0chn "
          "pDockQ pDockQ2 LIS n0res n0chn n0dom d0res d0chn d0dom nres1 nres2 dist1 dist2 Model\n")
def row(c1, c2, ipsae):
    f = [c1, c2, "10", "10", "asym", f"{ipsae:.6f}"] + ["0.0"] * 7 + ["50", "199", "0"] + ["0.0"] * 8 + ["m"]
    return " ".join(f) + "\n"

fails = []
def check(label, cond, detail=""):
    print(f"{'PASS' if cond else 'FAIL'}  {label}{(' -- ' + detail) if detail and not cond else ''}")
    if not cond: fails.append(label)

# ---- FAULT 2: a single-direction file must not yield a score --------------
with tempfile.TemporaryDirectory() as td:
    one = Path(td) / "one_direction_10_10.txt"
    one.write_text(HEADER + row("A", "B", 0.7654))
    two = Path(td) / "two_direction_10_10.txt"
    two.write_text(HEADER + row("A", "B", 0.7654) + row("B", "A", 0.8123))

    check("instrument_v2 refuses a one-direction file",
          SEEDAGG.read_cached(one)[0] is None, f"returned {SEEDAGG.read_cached(one)[0]!r}")
    check("master_rank refuses a one-direction file",
          CANON.ipsae_min(one) is None, f"returned {CANON.ipsae_min(one)!r}")
    check("instrument_v2 still scores a two-direction file",
          SEEDAGG.read_cached(two)[0] == 0.7654, f"returned {SEEDAGG.read_cached(two)[0]!r}")
    check("master_rank still scores a two-direction file",
          CANON.ipsae_min(two) == 0.7654, f"returned {CANON.ipsae_min(two)!r}")
    check("both take the MIN, not the first row or the max",
          SEEDAGG.read_cached(two)[0] == CANON.ipsae_min(two) == min(0.7654, 0.8123))

    # a 3-chain file: undefined without naming the pair, must refuse
    three = Path(td) / "three_pair_10_10.txt"
    three.write_text(HEADER + row("A", "B", 0.81) + row("B", "A", 0.82)
                            + row("A", "C", 0.22) + row("C", "A", 0.23)
                            + row("B", "C", 0.21) + row("C", "B", 0.24))
    check("instrument_v2 refuses a 3-pair file", SEEDAGG.read_cached(three)[0] is None)
    check("master_rank refuses a 3-pair file", CANON.ipsae_min(three) is None,
          f"returned {CANON.ipsae_min(three)!r} (min across INTERFACES, not directions)")

# ---- FAULT 1: a failed run must not return the previous run's number -----
src = HERE / "cases" / "01_barnase_barstar_positive_regression"
with tempfile.TemporaryDirectory() as td:
    d = Path(td) / "case"; shutil.copytree(src, d)
    cif = next(d.glob("*.cif")); pae = next(d.glob("*_ipsae.json"))
    good = LIVE.score(pae, cif)
    check("live path scores the intact case", good is not None and abs(good[0] - 0.888672) < 1e-6,
          f"got {good!r}")
    stale = cif.with_name(f"{cif.stem}_10_10.txt")
    check("the scored output exists on disk (the stale candidate)", stale.exists())
    before = stale.read_text()
    # corrupt the PAE so the reference cannot score, leaving the old output in place
    pae.write_text(json.dumps({"pae": "not-a-matrix"}))
    stale.write_text(before)                      # the previous run's answer, still there
    after = LIVE.score(pae, cif)
    check("live path refuses after the reference fails, despite a stale output file",
          after is None, f"returned {after[0] if after else after!r} from the stale file")

print()
if fails:
    print(f"{len(fails)} regression(s) FAILED: {', '.join(fails)}")
    sys.exit(1)
print("all fail-closed regressions pass: neither fault can return a number.")
