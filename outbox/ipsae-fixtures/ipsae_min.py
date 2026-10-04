#!/usr/bin/env python3
"""Compute ipSAE_min for predicted complexes, max over seeds.

ipsae.py reports per-direction rows plus a `max` row. The metric that predicts
binding (Overath et al.; Anthropic's protocol) is ipSAE_MIN -- the minimum over
both alignment directions. Taking the tool's `max` row is the easy mistake.

Usage:
  ipsae_min.py <pae_json_or_npz> <structure.cif> [pae_cutoff] [dist_cutoff]
  ipsae_min.py --dir <dir>     # every *_ipsae.json next to its *.cif
"""
import subprocess, sys, glob, os, statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
IPSAE = ROOT / "ipsae" / "ipsae.py"
PY = ROOT / ".venv" / "bin" / "python"


def score(pae_file: Path, struct: Path, pae_cut=10, dist_cut=10):
    """Run ipsae.py and return (ipsae_min, per_direction dict) or None."""
    # BUGFIX: resolve to absolute paths -- ipsae.py is run with cwd=struct.parent,
    # so relative paths from the caller's cwd would not resolve and it exits silently.
    pae_file, struct = pae_file.resolve(), struct.resolve()
    r = subprocess.run([str(PY), str(IPSAE), str(pae_file), str(struct),
                        str(pae_cut), str(dist_cut)],
                       capture_output=True, text=True, cwd=struct.parent)
    out = struct.with_name(f"{struct.stem}_{pae_cut}_{dist_cut}.txt")
    if not out.exists() and r.stderr:
        print(f"  ipsae stderr: {r.stderr.strip()[:300]}")
    if not out.exists():
        return None
    rows = [l.split() for l in out.read_text().splitlines()
            if l.strip() and not l.startswith("Chn1")]
    asym = {f"{r[0]}->{r[1]}": float(r[5]) for r in rows if len(r) > 5 and r[4] == "asym"}
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
    return min(next(iter(pairs.values()))), asym


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
