#!/usr/bin/env python3
"""pKa-perturbation sensitivity: does the submission's order survive PROPKA's own error?

The reviewer's direction on the reordered submission, 2026-10-05, restated: wherever a
candidate's selection flips, the deletion estimate has to be checked both against
alternatives prepared the same way (apo, relaxed) AND AGAINST pKa VALUES PUSHED AROUND
WITHIN THEIR PLAUSIBLE ERROR, with the per-candidate results kept rather than collapsed.
The output is a SENSITIVITY ANALYSIS and must not be dressed up as a measured confidence
interval; a ranking that shifts appreciably gets provisional tiers.

The three-basis comparison in METHODS 11.7 answered a different question -- how much the
COMPOSITION RULE moves the answer. This one asks how much PROPKA'S OWN UNCERTAINTY moves
it, which is the half of the ask addressed here. The apo/relaxed half needs new folding
and is recorded as unresolved.

METHOD. Every site's pKa_free and pKa_bound is already stored per pose in
ph_sensitivity.json. Each draw perturbs BOTH by independent Gaussian noise, recomputes
each site's linkage ratio, re-takes the histidine-only product (the ranking basis), takes
the median over that design's poses, and re-ranks all designs. Over many draws this gives
the distribution of each design's RANK.

SIGMA is PROPKA's own reported accuracy, not a tuned parameter. PROPKA 3 reports RMSD
~0.8 pKa units overall, worse for buried residues; 0.8 is used as the default and 0.4 and
1.2 are reported alongside so the conclusion can be read off any of them.

WHAT IS AND IS NOT REPORTED. Rank stability and how often pairs swap. NOT a confidence
interval on any ratio: the draws describe propagated pKa uncertainty only, and say nothing
about the fixed-conformation free leg, the independent-titration assumption behind
multiplying site ratios, or whether PROPKA is biased rather than noisy on these sites. A
consistent bias would move every design together and this analysis would not see it.

    ph_pka_perturbation.py [--sigma 0.8] [--draws 500]
"""
import csv, json, math, random, statistics as st, sys

SENS = 'analysis/01-egfr/ph_sensitivity.json'
CSV = 'submissions/01-egfr.csv'
OUT = 'analysis/01-egfr/ph_pka_perturbation.json'
PH_LO, PH_HI = 6.5, 7.4


def link(free, bound):
    """K(6.5)/K(7.4) for one proton-binding site. Same form as bin/ph_gate_all.link."""
    def K(ph):
        return (1 + 10 ** (bound - ph)) / (1 + 10 ** (free - ph))
    return K(PH_LO) / K(PH_HI)


def main():
    sigma = 0.8
    draws = 500
    if '--sigma' in sys.argv:
        sigma = float(sys.argv[sys.argv.index('--sigma') + 1])
    if '--draws' in sys.argv:
        draws = int(sys.argv[sys.argv.index('--draws') + 1])
    rng = random.Random(20261005)

    ship = {r['sequence'].strip().upper(): r['name'] for r in csv.DictReader(open(CSV))}
    sens = json.load(open(SENS))
    designs = {}
    for v in sens.values():
        q = v.get('seq', '').strip().upper()
        if q not in ship or 'sites_all_poses' not in v:
            continue
        # keep only the HISTIDINES -- the ranking basis
        poses = []
        for pose, sites in v['sites_all_poses'].items():
            his = [(s['free'], s['bound']) for s in sites.values() if s['resname'] == 'HIS']
            if his:
                poses.append(his)
        if poses:
            designs[ship[q]] = poses
    if not designs:
        sys.exit("no designs with per-site pKa data; re-run bin/ph_sensitivity_multisite.py")
    print(f"{len(designs)} submitted designs with per-site pKa data; "
          f"sigma = {sigma} pKa units, {draws} draws\n")

    # the unperturbed baseline, recomputed from the same stored pKa values
    def product(poses, jitter):
        meds = []
        for his in poses:
            p = 1.0
            for f, b in his:
                if jitter:
                    f = f + rng.gauss(0, sigma)
                    b = b + rng.gauss(0, sigma)
                p *= link(f, b)
            meds.append(p)
        return st.median(meds)

    base = {n: product(p, False) for n, p in designs.items()}
    base_order = [n for n, _ in sorted(base.items(), key=lambda kv: -kv[1])]
    base_rank = {n: i for i, n in enumerate(base_order, 1)}

    rank_draws = {n: [] for n in designs}
    top3_hits = {n: 0 for n in designs}
    for _ in range(draws):
        vals = {n: product(p, True) for n, p in designs.items()}
        order = [n for n, _ in sorted(vals.items(), key=lambda kv: -kv[1])]
        for i, n in enumerate(order, 1):
            rank_draws[n].append(i)
            if i <= 3:
                top3_hits[n] += 1

    print(f"{'design':<44}{'base':>6}{'rank':>7}{'p5-p95':>11}{'top3%':>7}")
    rows = []
    for n in base_order:
        rs = sorted(rank_draws[n])
        p5, p95 = rs[int(0.05 * len(rs))], rs[int(0.95 * len(rs)) - 1]
        med = st.median(rs)
        rows.append(dict(name=n, base_ratio=round(base[n], 4), base_rank=base_rank[n],
                         median_rank=med, p5_rank=p5, p95_rank=p95,
                         top3_fraction=round(top3_hits[n] / draws, 3)))
        print(f"{n[:43]:<44}{base_rank[n]:>6}{med:>7.0f}{f'{p5}-{p95}':>11}"
              f"{100*top3_hits[n]/draws:>7.0f}")

    # how often does the realised order match the baseline order?
    exact = sum(1 for n in designs if st.median(rank_draws[n]) == base_rank[n])
    moved = [r for r in rows if r['p95_rank'] - r['p5_rank'] >= 5]
    print(f"\n{exact} of {len(designs)} designs keep their baseline rank as the perturbed median.")
    print(f"{len(moved)} design(s) span 5 or more ranks across the central 90% of draws"
          + (f": {', '.join(r['name'][:28] for r in moved)}" if moved else ""))
    json.dump(dict(sigma=sigma, draws=draws, designs=rows), open(OUT, 'w'), indent=1)
    print(f"wrote {OUT}")


if __name__ == '__main__':
    main()
