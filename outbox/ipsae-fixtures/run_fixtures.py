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
IPSAE = REPO / "ipsae" / "ipsae.py"
PY_EXE = REPO / ".venv" / "bin" / "python"
if not PY_EXE.exists(): PY_EXE = sys.executable
PAE_CUT, DIST_CUT = 10, 10


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
                         "naming the binder:target pair", "asym": asym}
    import gemmi
    st = gemmi.read_structure(str(cif)); st.setup_entities()
    return {"structure": cif.name,
            "chains": {c.name: len(c) for c in st[0]},
            "asym": asym,
            "n_interface_residues": nres,
            "n_chainpair_residues": ntot,
            "ipsae_min": min(next(iter(pairs.values())))}


def main():
    res = {}
    for d in sorted(p for p in (HERE / "cases").iterdir() if p.is_dir()):
        res[d.name] = score_case(d)
        r = res[d.name]
        print(f"\n=== {d.name}")
        if "error" in r:
            print(f"    ERROR: {r['error']}"); continue
        print(f"    {r['structure']}")
        print(f"    chains {r['chains']}")
        for k, v in r["asym"].items():
            ni = r["n_interface_residues"].get(k, "?")
            nt = r["n_chainpair_residues"].get(k, "?")
            note = "   <-- ZERO interface residues: no confident interface, not a crash" if ni == "0" else ""
            print(f"    {k:>8}  ipSAE {v:.6f}   n0res (interface) {ni:>5}  of n0chn {nt}{note}")
        print(f"    ipSAE_min (min over DIRECTIONS) = {r['ipsae_min']:.6f}")
    exp = HERE / "expected.json"
    if "--check" in sys.argv:
        want = json.loads(exp.read_text()); bad = 0
        for k, v in want.items():
            got = res.get(k, {}).get("ipsae_min")
            ok = got is not None and abs(got - v["ipsae_min"]) < 1e-6
            print(f"{'OK  ' if ok else 'FAIL'} {k}: expected {v['ipsae_min']:.6f} got {got}")
            bad += (not ok)
        sys.exit(1 if bad else 0)
    exp.write_text(json.dumps(res, indent=1))
    print(f"\nwrote {exp}")


if __name__ == "__main__":
    main()
