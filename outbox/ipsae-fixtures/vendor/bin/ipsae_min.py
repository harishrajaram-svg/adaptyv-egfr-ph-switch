#!/usr/bin/env python3
"""Compute ipSAE_min for predicted complexes, max over seeds.

ipsae.py reports per-direction rows plus a `max` row. The metric that predicts
binding (Overath et al.; Anthropic's protocol) is ipSAE_MIN -- the minimum over
both alignment directions. Taking the tool's `max` row is the easy mistake.

Usage:
  ipsae_min.py <pae_json_or_npz> <structure.cif> [pae_cutoff] [dist_cutoff]
  ipsae_min.py --dir <dir>     # every *_ipsae.json next to its *.cif
"""
import subprocess, sys, glob, os, math, statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# THE REFERENCE MUST BE FINDABLE IN A CLONE.
#
# `/ipsae/` is gitignored, so the working-tree copy at ROOT/ipsae/ipsae.py is NOT published.
# The published copy is the pinned, sha256-guarded one vendored with the fixtures. Until
# 2026-10-05 this looked only at ROOT and therefore failed closed in every clone -- so the
# production scorer returned None for all 12 fixture cases for any reader who checked out
# the repository, while passing locally. That is the same class of fault the external
# reviewer reported as a 404 on the vendored reference, one level deeper: the file was
# published and the code still could not find it.
#
# Order matters: the working tree wins when present, so local runs are unchanged, and the
# vendored copy is the fallback rather than the default.
_V = ROOT / "outbox" / "ipsae-fixtures" / "vendor" / "ipsae" / "ipsae.py"
IPSAE = ROOT / "ipsae" / "ipsae.py"
if not IPSAE.exists() and _V.exists():
    IPSAE = _V
# A CLONE HAS NO .venv. This was hardcoded, so the production scorer -- the module the
# external reviewer was asked to check, and the first thing any reader runs -- died with
# FileNotFoundError on '<clone>/.venv/bin/python' in a fresh checkout. Prefer the repo venv
# when it is there, otherwise use the interpreter already running, which is what a reader
# invoking this script actually has.
_PY = ROOT / ".venv" / "bin" / "python"
PY = _PY if _PY.exists() else Path(sys.executable)


def score(pae_file: Path, struct: Path, pae_cut=10, dist_cut=10):
    """Run ipsae.py and return (ipsae_min, per_direction dict) or None.

    Fails closed: a score comes back only when ALL of these hold.
      1. the scoring process exits 0;
      2. the output file was written by THIS invocation -- any earlier output is
         deleted first, so a crashed run cannot resurrect a previous score;
      3. exactly one inter-chain pair is present (see the pair note below);
      4. that pair carries BOTH reciprocal directions, each finite.
    (4) matters because ipSAE_min is the min over two alignment directions. A
    file holding only one direction used to pass, reporting that single
    direction as if it were a minimum.
    """
    # BUGFIX: resolve to absolute paths -- ipsae.py is run with cwd=struct.parent,
    # so relative paths from the caller's cwd would not resolve and it exits silently.
    pae_file, struct = pae_file.resolve(), struct.resolve()
    out = struct.with_name(f"{struct.stem}_{pae_cut}_{dist_cut}.txt")
    if out.exists():
        out.unlink()                      # (2) invalidate before running
    r = subprocess.run([str(PY), str(IPSAE), str(pae_file), str(struct),
                        str(pae_cut), str(dist_cut)],
                       capture_output=True, text=True, cwd=struct.parent)
    if r.returncode != 0:                 # (1) non-zero exit is a failure, period
        msg = (r.stderr or r.stdout or "").strip()[:300]
        print(f"  ipsae exit {r.returncode}: {msg}")
        return None
    if not out.exists():
        if r.stderr:
            print(f"  ipsae stderr: {r.stderr.strip()[:300]}")
        return None
    rows = [l.split() for l in out.read_text().splitlines()
            if l.strip() and not l.startswith("Chn1")]
    asym = {}
    for row in rows:
        if len(row) > 5 and row[4] == "asym":
            try:
                v = float(row[5])
            except ValueError:
                continue
            if math.isfinite(v):
                asym[f"{row[0]}->{row[1]}"] = v
    if not asym:
        return None
    # ipSAE_min is the min over the TWO ALIGNMENT DIRECTIONS of ONE interface.
    # min() over the flat dict instead collapses across chain PAIRS, so on a
    # 3-chain complex (e.g. the g532 ladder: A=EGFR, B=VH, C=VL) the intra-binder
    # VH:VL interface enters the pool and a badly packed Fv can be reported as the
    # binder's binding score. With two chains there is one pair and this is
    # bit-identical to min(asym.values()) -- verified on 536 cached outputs.
    pairs = {}
    for k, v in asym.items():
        c1, c2 = k.split("->")
        pairs.setdefault(frozenset((c1, c2)), []).append(v)
    if len(pairs) > 1:
        pretty = ", ".join(f"{':'.join(sorted(p))}={min(v):.4f}" for p, v in pairs.items())
        raise SystemExit(
            f"{struct}: {len(pairs)} inter-chain pairs ({pretty}).\n"
            "  ipSAE_min is undefined without naming the binder:target pair; refusing\n"
            "  to let min() collapse across interfaces. Score the pair explicitly.")
    pair, vals = next(iter(pairs.items()))
    if len(pair) != 2 or len(vals) != 2:  # (4) exactly two reciprocal directions
        c = ":".join(sorted(pair))
        print(f"  {struct.name}: pair {c} has {len(vals)} finite asym direction(s), "
              f"need 2 -- ipSAE_min is a min over two directions, not a single one")
        return None
    return min(vals), asym


def main():
    args = sys.argv[1:]
    if args and args[0] == "--dir":
        pairs = []
        for j in sorted(glob.glob(os.path.join(args[1], "**", "*_ipsae.json"), recursive=True)):
            cif = Path(j.replace("_ipsae.json", ".cif"))
            if cif.exists():
                pairs.append((Path(j), cif))
        if not pairs:
            sys.exit(f"no *_ipsae.json + .cif pairs under {args[1]}")
    else:
        pairs = [(Path(args[0]), Path(args[1]))]

    results = {}
    for pae, cif in pairs:
        r = score(pae, cif)
        if r is None:
            print(f"FAILED  {cif.name}  (no ipSAE output -- is the PAE matrix present?)")
            continue
        m, asym = r
        # strip trailing _seed<N>_sample_<N> so all seeds of one design group together
        import re as _re
        design = _re.sub(r"(_seed\d+)?_sample_\d+$", "", cif.stem)
        results.setdefault(design, []).append(m)
        print(f"{cif.name}: ipSAE_min={m:.4f}  ({', '.join(f'{k}={v:.3f}' for k, v in asym.items())})")

    if results:
        print("\n=== per design, MAX over seeds (the protocol aggregation) ===")
        for d, vals in sorted(results.items(), key=lambda kv: -max(kv[1])):
            spread = max(vals) - min(vals)
            extra = (f", spread={spread:.4f}, min={min(vals):.4f}" if len(vals) > 1 else "")
            print(f"  {d}: ipSAE_min={max(vals):.4f}  (n={len(vals)}{extra})")


if __name__ == "__main__":
    main()
