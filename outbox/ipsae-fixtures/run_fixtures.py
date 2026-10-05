#!/usr/bin/env python3
"""Reproduce every ipSAE_min fixture. No arguments.

    python3 run_fixtures.py            # score every case, print both directions + the min
    python3 run_fixtures.py --check    # compare against expected.json, exit 1 on mismatch

Scores each case with the PINNED Dunbrack reference (../../ipsae/ipsae.py, commit
6174cf9e71cb1bd660cc805856a18c4871a6dec3) and reports BOTH asymmetric directions
separately before taking the minimum, so the min-over-DIRECTIONS step is visible and
auditable rather than buried.

NOTE ON A DISTINCTION PK RAISED: the reference's own `max` row is a maximum across chain
DIRECTIONS of one interface. That is a different question from our max-versus-median across
SEEDS. This script only concerns the former. Seed aggregation happens in
bin/instrument_v2.py and is reported separately.
"""
import json, os, subprocess, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
# Resolve the pinned Dunbrack reference. The bundle VENDORS it under vendor/ipsae/ (MIT,
# commit 6174cf9e71cb1bd660cc805856a18c4871a6dec3) so this runs on a machine that has only
# this directory. The repo checkout is preferred when present so a reviewer can confirm the
# two copies are identical; `ipsae/` is gitignored in our repo, which is why the vendored
# copy exists at all -- without it every case errored out on anyone else's machine.
IPSAE = HERE / "vendor" / "ipsae" / "ipsae.py"
if (REPO / "ipsae" / "ipsae.py").exists():
    IPSAE = REPO / "ipsae" / "ipsae.py"
if not IPSAE.exists():
    sys.exit(f"reference implementation not found at {IPSAE}")
# A CLONE HAS NO .venv -- see the same fix in bin/gate_sweep.py. Hardcoding it meant the
# fixture bundle, whose entire purpose is to let a reviewer reproduce the scorer, could not
# run in a fresh checkout. PK reported exactly this class of problem for the vendored
# reference (a 404); this was the same failure one level up.
import sys as _sys
_PV = REPO / ".venv" / "bin" / "python"
PY_EXE = _PV if _PV.exists() else Path(_sys.executable)
if not PY_EXE.exists(): PY_EXE = sys.executable
PAE_CUT, DIST_CUT = 10, 10


# ---------------------------------------------------------------------------
# THE PRODUCTION PARSERS, loaded from bin/ and exercised on every case.
#
# PK, 2026-10-05: "The runner invokes the reference and does its own parsing --
# it does NOT exercise bin/ipsae_min.py." That was correct, and it was the whole
# weakness of this bundle: it proved the REFERENCE reproduces, not that OUR code
# reads it correctly. Three separate copies of the parse logic existed (bin/
# ipsae_min.py for live scoring, bin/instrument_v2.py for seed aggregation, and
# bin/master_rank.py for the canonical file the submission is built from), and
# the two faults PK found lived in the copies, not in this runner's.
#
# Every case is now scored FOUR ways -- the reference, plus all three production
# parsers -- and --check requires all four to agree. A divergence between the
# copies is now a test failure instead of an invisible inconsistency.
import importlib.util

def _load(name):
    """Load a production parser. The repo checkout is preferred so a reviewer sees the
    live code; vendor/bin/ is the self-contained fallback for the standalone bundle.

    This HARD-FAILS when neither exists. Returning None would make --check quietly
    skip the production comparison and still print "11/11 cases reproduce" -- the
    exact fail-open shape of the two faults this bundle was rebuilt to catch."""
    path = REPO / "bin" / f"{name}.py"
    if not path.exists():
        path = HERE / "vendor" / "bin" / f"{name}.py"
    if not path.exists():
        sys.exit(f"production parser {name}.py not found in {REPO/'bin'} or "
                 f"{HERE/'vendor'/'bin'} -- refusing to run a check that would skip it")
    spec = importlib.util.spec_from_file_location(f"_prod_{name}", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)          # safe: all three are __main__-guarded
    return mod

PROD_LIVE     = _load("ipsae_min")        # .score(pae, cif)  -- runs the reference itself
PROD_SEEDAGG  = _load("instrument_v2")    # .read_cached(txt)
PROD_CANONICAL= _load("master_rank")      # .ipsae_min(txt)   -- feeds master_rank.json


def production_scores(cif: Path, pae: Path, out: Path):
    """Score one case through each production code path. None = that path refused."""
    got = {}
    if True:
        try:
            r = PROD_LIVE.score(pae, cif, PAE_CUT, DIST_CUT)   # regenerates `out`
            got["bin/ipsae_min.py"] = None if r is None else r[0]
        except SystemExit as e:
            got["bin/ipsae_min.py"] = f"refused: {str(e).splitlines()[0][:80]}"
    if True:
        got["bin/instrument_v2.py"] = (PROD_SEEDAGG.read_cached(out)[0]
                                       if out.exists() else None)
    if True:
        got["bin/master_rank.py"] = (PROD_CANONICAL.ipsae_min(out)
                                     if out.exists() else None)
    return got


def score_case(d: Path):
    cif = next(iter(sorted(d.glob("*.cif"))), None)
    pae = next(iter(sorted(d.glob("*_ipsae.json"))), None)
    if cif is None or pae is None:
        return {"error": "missing cif or pae"}
    out = cif.with_name(f"{cif.stem}_{PAE_CUT}_{DIST_CUT}.txt")
    if out.exists(): out.unlink()
    r = subprocess.run([str(PY_EXE), str(IPSAE), str(pae.resolve()), str(cif.resolve()),
                        str(PAE_CUT), str(DIST_CUT)],
                       capture_output=True, text=True, cwd=str(d))
    if not out.exists():
        return {"error": "ipsae.py produced no output", "stderr": r.stderr.strip()[:300]}
    rows = [l.split() for l in out.read_text().splitlines()
            if l.strip() and not l.startswith("Chn1")]
    asym = {f"{x[0]}->{x[1]}": float(x[5]) for x in rows if len(x) > 5 and x[4] == "asym"}
    # column 13 is n0res -- the count of residue pairs PASSING the PAE/dist interface filter.
    # Column 14 is n0chn, the total chain-pair residue count, which is NOT the interface size.
    # An earlier version of this script printed column 14 and mislabelled it "interface
    # residues"; that made a zero-interface case look like a full interface scoring zero.
    nres = {f"{x[0]}->{x[1]}": x[13] for x in rows if len(x) > 14 and x[4] == "asym"}
    ntot = {f"{x[0]}->{x[1]}": x[14] for x in rows if len(x) > 14 and x[4] == "asym"}
    # min over the two ALIGNMENT DIRECTIONS of ONE interface -- never across chain PAIRS
    pairs = {}
    for k, v in asym.items():
        a, b = k.split("->"); pairs.setdefault(frozenset((a, b)), []).append(v)
    if len(pairs) != 1:
        return {"error": f"{len(pairs)} inter-chain pairs; ipSAE_min is undefined without "
                         "naming the binder:target pair", "asym": asym,
                "production": production_scores(cif, pae, out)}
    # gemmi is used only to report chain sizes in the printout. It is NOT on the
    # scoring path, so a reviewer without it still reproduces every number.
    try:
        import gemmi
        st = gemmi.read_structure(str(cif)); st.setup_entities()
        chains = {c.name: len(c) for c in st[0]}
    except ImportError:
        chains = {"(gemmi not installed; chain sizes omitted)": 0}
    return {"structure": cif.name,
            "chains": chains,
            "asym": asym,
            "n_interface_residues": nres,
            "n_chainpair_residues": ntot,
            "ipsae_min": min(next(iter(pairs.values()))),
            "production": production_scores(cif, pae, out)}


def check_vendor_drift():
    """vendor/bin/ must be byte-identical to bin/ when both are present.

    A vendored copy that drifts is how outbox/ipsae-fixtures/ipsae_min.py came to be a
    stale fork of bin/ipsae_min.py: the bundle tested one file while the submission was
    built by another. Fail loudly rather than test the wrong code."""
    import hashlib
    drift = []
    for name in ("ipsae_min", "instrument_v2", "master_rank"):
        live, vend = REPO / "bin" / f"{name}.py", HERE / "vendor" / "bin" / f"{name}.py"
        if live.exists() and vend.exists():
            h = lambda f: hashlib.sha256(f.read_bytes()).hexdigest()
            if h(live) != h(vend):
                drift.append(f"  {name}.py: bin/={h(live)[:12]} vendor/bin/={h(vend)[:12]}")
    if drift:
        sys.exit("vendor/bin has drifted from bin/:\n" + "\n".join(drift) +
                 "\n  Re-copy before trusting this bundle:\n"
                 "    cp bin/{ipsae_min,instrument_v2,master_rank}.py "
                 "outbox/ipsae-fixtures/vendor/bin/")


def main():
    check_vendor_drift()
    res = {}
    for d in sorted(p for p in (HERE / "cases").iterdir() if p.is_dir()):
        res[d.name] = score_case(d)
        r = res[d.name]
        print(f"\n=== {d.name}")
        if "error" in r:
            print(f"    ERROR: {r['error']}")
            for path, val in (r.get("production") or {}).items():
                shown = f"{val:.6f}" if isinstance(val, float) else repr(val)
                mark = "refuses too" if val is None or isinstance(val, str) else "RETURNED A NUMBER -- fail-open"
                print(f"      via {path:<24} {shown:>12}   {mark}")
            continue
        print(f"    {r['structure']}")
        print(f"    chains {r['chains']}")
        for k, v in r["asym"].items():
            ni = r["n_interface_residues"].get(k, "?")
            nt = r["n_chainpair_residues"].get(k, "?")
            note = "   <-- ZERO interface residues: no confident interface, not a crash" if ni == "0" else ""
            print(f"    {k:>8}  ipSAE {v:.6f}   n0res (interface) {ni:>5}  of n0chn {nt}{note}")
        print(f"    ipSAE_min (min over DIRECTIONS) = {r['ipsae_min']:.6f}")
        for path, val in (r.get("production") or {}).items():
            mark = "agrees" if isinstance(val, float) and abs(val - r["ipsae_min"]) < 1e-9 else "DIVERGES"
            shown = f"{val:.6f}" if isinstance(val, float) else repr(val)
            print(f"      via {path:<24} {shown:>12}   {mark}")
    exp = HERE / "expected.json"
    if "--check" in sys.argv:
        want = json.loads(exp.read_text()); bad = 0
        print()
        for k, v in want.items():
            got = res.get(k, {})
            # A case whose EXPECTED result is a refusal is checked on the refusal, not on a
            # score. Case 10 is exactly that: a 3-chain complex where ipSAE_min is undefined
            # until the binder:target pair is named, and the correct behaviour is to refuse.
            # The earlier version of this loop assumed every case yields a number and crashed
            # with KeyError on the refusal case -- i.e. the check could not express "the right
            # answer here is an error", which is the one behaviour PK asked us to demonstrate
            # ("a failed run must not silently become a valid score of zero").
            if "error" in v:
                ok = "error" in got and got["error"] == v["error"]
                print(f"{'OK  ' if ok else 'FAIL'} {k}: expected REFUSAL -> "
                      f"{'refused as expected' if ok else got.get('ipsae_min', got.get('error', 'no result'))}")
            else:
                g = got.get("ipsae_min")
                ok = g is not None and abs(g - v["ipsae_min"]) < 1e-6
                print(f"{'OK  ' if ok else 'FAIL'} {k}: expected {v['ipsae_min']:.6f} got {g}")
                # BOTH DIRECTIONAL VALUES, not just the minimum.
                #
                # Reviewer, 2026-10-05: "Its check compares the final minimum, not both
                # directional values." Correct -- the directions were printed and never
                # asserted. Two different (A->B, B->A) pairs can share a minimum, so a
                # check on the min alone passes while the asymmetry underneath it has
                # changed. That asymmetry is the whole quantity this bundle exists to
                # document: case 03 and case 12 differ from their partners ONLY in the
                # direction that is not the minimum.
                for d, want_v in (v.get("asym") or {}).items():
                    gv = (got.get("asym") or {}).get(d)
                    if gv is None:
                        print(f"     FAIL {k}: direction {d} missing from the result")
                        ok = False
                    elif abs(gv - want_v) > 1e-6:
                        print(f"     FAIL {k}: direction {d} expected {want_v:.6f} got {gv:.6f}")
                        ok = False
                extra = set((got.get("asym") or {})) - set((v.get("asym") or {}))
                if extra:
                    print(f"     FAIL {k}: unexpected direction(s) {sorted(extra)}")
                    ok = False
            # EVERY production parser must land on the reference answer. A case that
            # reproduces against the reference while our own code reads it differently
            # is not a passing case -- that gap is what this bundle exists to close.
            want_num = v.get("ipsae_min")
            for path, val in (got.get("production") or {}).items():
                if want_num is None:                  # refusal expected
                    pok = val is None or isinstance(val, str)
                    detail = "refused" if pok else f"returned {val!r} (fail-open)"
                else:
                    pok = isinstance(val, float) and abs(val - want_num) < 1e-6
                    detail = (f"{val:.6f}" if isinstance(val, float) else repr(val))
                if not pok:
                    print(f"     FAIL {k} via {path}: {detail}")
                    ok = False
            bad += (not ok)
        print(f"\n{len(want) - bad}/{len(want)} cases reproduce.")
        sys.exit(1 if bad else 0)
    if "--freeze" not in sys.argv and exp.exists():
        # expected.json is the frozen reference. Overwriting it on every plain run means a
        # regression silently becomes the new expectation -- the same class of error as a
        # threshold chosen after seeing the data. Writing it now requires --freeze.
        print(f"\n{exp.name} left unchanged. Re-freeze deliberately with --freeze.")
        return
    exp.write_text(json.dumps(res, indent=1))
    print(f"\nwrote {exp}")


if __name__ == "__main__":
    main()
