#!/usr/bin/env python3
"""Combine per-arm ipSAE_min into the protocol's ensemble ranking.

Anthropic's instrument: score ipSAE_min on each co-folding arm, take the MAX
over seeds, z-score WITHIN target across designs, then average the arms.
That 3-arm z-ensemble reached 0.66 macro-AP vs 0.55 for AlphaFold3 alone.

Two rules this encodes, both easy to get wrong:
  * Raw scores are NOT comparable across targets. z-scores are transductive --
    they only mean something within one target's pool of designs. Never carry a
    z-score between targets, waves or campaigns.
  * Ranking needs a POOL. A single design has no z-score; with n<3 this reports
    raw means and says so.

Input: JSON {"<arm>": {"<design>": ipsae_min, ...}, ...}
Usage: ensemble_rank.py scores.json [--weights arm=4,arm2=1]
"""
import json, sys, statistics as st


def zscores(vals: dict):
    xs = list(vals.values())
    if len(xs) < 3:
        return None
    mu = st.mean(xs)
    sd = st.pstdev(xs)
    if sd == 0:
        return {k: 0.0 for k in vals}
    return {k: (v - mu) / sd for k, v in vals.items()}


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    arms = json.load(open(sys.argv[1]))

    weights = {a: 1.0 for a in arms}
    for a in sys.argv[2:]:
        if a.startswith("--weights"):
            for kv in a.split("=", 1)[1].split(","):
                k, v = kv.split("=")
                weights[k] = float(v)

    designs = sorted({d for arm in arms.values() for d in arm})
    n = len(designs)
    zs = {a: zscores(arms[a]) for a in arms}
    transductive = all(z is not None for z in zs.values())

    if not transductive:
        print(f"WARNING: only {n} design(s). z-scoring needs a pool (n>=3).")
        print("Reporting the raw weighted mean instead -- NOT comparable to any other target.\n")

    rows = []
    for d in designs:
        parts, wsum = 0.0, 0.0
        for a in arms:
            if d not in arms[a]:
                continue
            v = zs[a][d] if transductive else arms[a][d]
            w = weights.get(a, 1.0)
            parts += w * v
            wsum += w
        rows.append((d, parts / wsum if wsum else float("nan"),
                     {a: arms[a].get(d) for a in arms}))

    rows.sort(key=lambda r: -r[1])
    label = "ensemble_z" if transductive else "weighted_raw"
    print(f"{'rank':<5}{'design':<34}{label:>14}   per-arm ipSAE_min")
    for i, (d, s, per) in enumerate(rows, 1):
        detail = "  ".join(f"{a}={v:.3f}" if v is not None else f"{a}=--"
                           for a, v in per.items())
        print(f"{i:<5}{d:<34}{s:>14.4f}   {detail}")

    if transductive:
        print(f"\nz-scored within this target across {n} designs. Do not compare to another target.")


if __name__ == "__main__":
    main()
