#!/usr/bin/env python3
"""Combined pH gate: EVERY titratable site, on BOTH partners, in one pass.

WHY THIS FILE EXISTS, AND WHY IT IS THE MOST DANGEROUS CODE IN THIS PROJECT
--------------------------------------------------------------------------
Two gates existed and neither could see a two-site design:

    bin/ph_gate_refolds.py / ph_gate_all.py   TARGET histidines only (chain = target)
    bin/ph_gate_mechA.py                      BINDER histidines only (chain = binder)

They were never run on the same molecule. A genuinely two-site design -- a carboxylate on
the binder shifting target H433, PLUS a binder histidine paired to a target carboxylate --
would come back from the first gate reporting H433 alone and from the second reporting the
binder histidine alone, with the linkage between them invisible in both.

This script composes them. Its entire purpose is therefore to produce a number ABOVE the
single-site ceiling, and METHODS section 7 records that we have already manufactured exactly
such a number once, by composing over the sites that helped:

    "`product` composed over helping sites only -- selection on the outcome; manufactured
     our only claim above the 7.94x ceiling"

rimA01/d3_rimA_50 was reported at 9.10x (H370 6.056 x H383 1.503) while its own H433 read
0.723. Its true all-site product is 6.581x. 1,547 of 1,584 designs were overstated, by up
to 4.77x. So this file is written to make that class of error structurally impossible, and
every guard below exists because of it. Read FIVE GUARDS before trusting any output.

THE PHYSICS
-----------
Thermodynamic linkage for one titratable site:

    K(pH) = (1 + 10^(pKa_bound - pH)) / (1 + 10^(pKa_free - pH))
    ratio = K(6.5) / K(7.4)

Independent sites multiply. The one-proton bound over this pH pair is 10^0.9 = 7.943x; for
a site with pKa_free = 6.22 the attainable range is 0.699x to 5.554x (both figures
independently confirmed by the reviewer, 2026-10-03).

THE FREE LEG, ON BOTH SIDES, BY DELETION IN PLACE
-------------------------------------------------
ph_gate_all computes the target's free pKa by deleting the BINDER from the same file,
holding every other coordinate fixed -- METHODS section 2 defends that at length, because
taking pKa_free from a separate apo structure turned 10 apparent switches into 1 when the
generator repacked the target per design.

This script applies the mirror image to the binder: its free pKa comes from the same file
with the TARGET deleted in place. That is symmetric, it needs no additional folding, and it
holds the binder's own conformation fixed. ph_gate_mechA instead required a separately
folded binder-only structure (--binder-dir), which introduces a conformational difference
the deletion method does not have.

Cost: 4 PROPKA runs per pose rather than 2. PROPKA is local CPU.

FIVE GUARDS
-----------
 1. ALL-SITES COMPOSITION, never a subset. Sites with ratio < 1 are sites where acid
    genuinely weakens binding and they belong in the product. `product_helpful_only` is
    computed and reported ONLY so the old bug's value stays visible for audit; it is never
    the headline and `--check` fails if anything reads it as one.
 2. CEILING CHECK. If k sites moved, the product cannot exceed 7.943^k. A product above
    that is arithmetically impossible and is reported as `IMPOSSIBLE`, not as a result.
 3. IMPLIED pKa_bound. Every ratio is inverted back to the pKa_bound it requires. A site
    demanding pKa_bound > 10.5 or < 2.0 for a histidine is flagged `implausible` -- a real
    histidine does not reach those, so the number is PROPKA noise or a packing artifact.
 4. COUNTER-CHARGE DISTANCE. A site with a large shift and no counter-charge within 6 A on
    the other chain is flagged `no_partner`. That is the H370 desolvation signature: a shift
    with no mechanism. Measured from the titratable atom to the nearest opposite-charge
    heavy atom ACROSS the interface.
 5. EPISTASIS RECONSTRUCTION. With --parent and --singles, the double mutant's product is
    checked against its parts. A product that cannot be reconstructed from the single
    mutants is epistasis or a bug, and the script says which it cannot distinguish.
    This is the check that would have caught the 9.10x.

USAGE
    ph_gate_multisite.py <complex.cif> [more...]            # score poses
    ph_gate_multisite.py --dir <run_dir> [--tsv out.tsv]    # all poses under a run
    ph_gate_multisite.py --selftest                         # no I/O, pure arithmetic
"""
import argparse, glob, json, math, os, re, subprocess, sys, tempfile
from pathlib import Path

# THE pH PAIR IS PER PROBLEM, AND GETTING IT FROM THE WRONG MODULE IS A SILENT WRONG ANSWER.
# Problem 1 is bind-at-6.5 / silent-at-7.4. Problem 2 is bind-at-7.4 / silent-at-6.0. Until
# 2026-10-08 this module held its own PH_LO/PH_HI while ph_gate_all.link() read ph_gate_all's,
# so setting the pair here for problem 2 would have left the physics on 6.5 with no error.
# link() now takes the pair, and EVERY call in this file goes through _link() so no call site
# can forget. The selftest greps this source to prove there is no bare link( left.
PROBLEMS = {
    # problem: (ph_lo, ph_hi, merit_direction, what the design must do)
    1: (6.5, 7.4, 'above_1', 'bind at 6.5, silent at 7.4 -- acid must TIGHTEN, ratio > 1'),
    2: (6.0, 7.4, 'below_1', 'bind at 7.4, silent at 6.0 -- acid must WEAKEN, ratio < 1'),
}
PROBLEM = 1                                    # default preserves every problem-1 result
PH_LO, PH_HI, MERIT, MERIT_TEXT = PROBLEMS[PROBLEM]
ONE_PROTON_BOUND = 10 ** (PH_HI - PH_LO)      # 7.943x at 6.5/7.4; 25.119x at 6.0/7.4


def set_problem(n):
    """Switch the pH pair and the ceiling together. They must never move apart."""
    global PROBLEM, PH_LO, PH_HI, MERIT, MERIT_TEXT, ONE_PROTON_BOUND
    if n not in PROBLEMS:
        sys.exit(f"REFUSE: unknown problem {n}; have {sorted(PROBLEMS)}")
    PROBLEM = n
    PH_LO, PH_HI, MERIT, MERIT_TEXT = PROBLEMS[n]
    ONE_PROTON_BOUND = 10 ** (PH_HI - PH_LO)


def orient_multi(model):
    """(target_chains, binder_chain, names) -- generalised from the two-chain case.

    ph_gate_all.orient() tries chain 'B' then 'A' and returns ONE target chain. That is right
    for problem 1, where the complex is target + binder. Problem 2's target is a HOMOTRIMER, so
    the complex is three target protomers plus a binder, and the single-target assumption breaks
    in a way that matters: the old free-target leg deleted every chain that was not the chosen
    target, which on a trimer removes the two SIBLING PROTOMERS along with the binder. H73 sits
    at the inter-protomer seam -- it is the residue BinderBench names as a hotspot for this
    target -- so deleting a sibling would hand it a free pKa from a monomer that does not exist
    in the assay. The free leg must remove the BINDER and nothing else.

    Identification is by fingerprint, never by chain order or count: every chain whose histidine
    tuple matches a registered family is a target protomer, and the one remaining chain is the
    binder. Refuses anything that is not N target chains plus exactly one binder.
    """
    chains = list(model)
    matched, names = [], None
    for c in chains:
        f = family(c)
        if f is not None:
            matched.append(c.name)
            names = f
    others = [c.name for c in chains if c.name not in matched]
    if names is None:
        return None, None, None, 'target family unrecognised'
    if len(others) != 1:
        return None, None, None, (f'{len(matched)} target chain(s) {matched} and '
                                  f'{len(others)} non-target {others}; need exactly one binder')
    return sorted(matched), others[0], names, None


def _link(free, bound):
    """The ONLY linkage path in this file. Always carries this module's pH pair."""
    return link(free, bound, PH_LO, PH_HI)


def merits(ratio):
    """Does this ratio point the way THIS problem needs? Direction is not cosmetic: a
    ratio of 5x is a success for problem 1 and a failure for problem 2."""
    if abs(ratio - 1.0) < NOISE:
        return 'noise'
    return 'helps' if ((ratio > 1.0) == (MERIT == 'above_1')) else 'wrong_direction'
NOISE = 0.05                                   # |ratio-1| below this is PROPKA noise
# A pKa_bound outside the window for that residue's chemistry is suspect. Acids get
# their own window: a carboxylate cannot titrate at a histidine's pKa.
PKA_PLAUSIBLE_BY_RES = {'HIS': (2.0, 10.5), 'ASP': (0.5, 9.0), 'GLU': (0.5, 9.0)}
PKA_PLAUSIBLE = PKA_PLAUSIBLE_BY_RES['HIS']    # kept: referenced by the selftest
PARTNER_CUT = 6.0                              # A, titratable atom -> nearest counter-charge

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ph_gate_all import CROP, D3, ECD, link, family, orient   # one source for the physics

ACID_O = {('ASP','OD1'),('ASP','OD2'),('GLU','OE1'),('GLU','OE2')}
BASE_N = {('HIS','ND1'),('HIS','NE2'),('LYS','NZ'),('ARG','NH1'),('ARG','NH2'),('ARG','NE')}


def implied_pka_bound(ratio, free):
    """Invert the linkage equation for pKa_bound. Returns None if the ratio is unreachable."""
    lo, hi = free - 8.0, free + 20.0
    if not (_link(free, lo) <= ratio <= _link(free, hi)):
        return None
    for _ in range(200):
        mid = (lo + hi) / 2
        if _link(free, mid) < ratio: lo = mid
        else: hi = mid
    return round((lo + hi) / 2, 2)


def propka_his_and_acids(st, wd, tag, chain):
    """pKa of every HIS/ASP/GLU on `chain`. Returns {(resname,resnum): pKa}."""
    st.write_pdb(os.path.join(wd, tag + '.pdb'))
    subprocess.run([sys.executable, '-m', 'propka', tag + '.pdb'], capture_output=True, cwd=wd)
    f = os.path.join(wd, tag + '.pka'); out = {}
    if os.path.exists(f):
        for ln in open(f):
            q = ln.split()
            if len(q) > 3 and q[0] in ('HIS', 'ASP', 'GLU') and q[2] == chain:
                try: out[(q[0], int(q[1]))] = float(q[3])
                except ValueError: pass
    return out


def cross_partner_distance(model, chain_of_site, resnum, resname, other_chain):
    """Nearest opposite-charge heavy atom on the OTHER chain. Guard 4."""
    import gemmi
    C = {c.name: c for c in model}
    if chain_of_site not in C or other_chain not in C: return None
    site = next((r for r in C[chain_of_site] if r.seqid.num == resnum and r.name == resname), None)
    if site is None: return None
    want_acid = resname in ('HIS', 'LYS', 'ARG')          # a base looks for acids
    mine = [a.pos for a in site if (resname, a.name) in (BASE_N if want_acid else ACID_O)]
    if not mine: return None
    theirs = [a.pos for r in C[other_chain] for a in r
              if (r.name, a.name) in (ACID_O if want_acid else BASE_N)]
    if not theirs: return None
    return round(min(p.dist(q) for p in mine for q in theirs), 2)


def score_pose(cif):
    """All titratable sites on both partners, one pose. Returns a dict or None."""
    import gemmi
    try:
        st = gemmi.read_structure(cif); st.setup_entities()
        st.remove_ligands_and_waters(); st.setup_entities()
        tgts, bnd, names, why = orient_multi(st[0])
        if why:
            return {'file': os.path.basename(cif), 'error': why}

        with tempfile.TemporaryDirectory() as wd:
            # BOUND legs: the intact complex, read per chain.
            bound_t = {}
            for c in tgts:
                for k, v in propka_his_and_acids(st, wd, f'cpx_t_{c}', c).items():
                    bound_t[(c,) + k] = v
            bound_b = propka_his_and_acids(st, wd, 'cpx_b', bnd)

            # FREE TARGET: delete the BINDER ONLY. On a homotrimer the sibling protomers stay,
            # because they are present in the assay and H73 titrates against them.
            a = gemmi.read_structure(cif); a.setup_entities()
            a.remove_ligands_and_waters(); a.setup_entities()
            for i in range(len(a[0]) - 1, -1, -1):
                if a[0][i].name == bnd: del a[0][i]
            a.setup_entities()
            assert len([c.name for c in a[0]]) == len(tgts), 'free-target leg lost a protomer'
            free_t = {}
            for c in tgts:
                for k, v in propka_his_and_acids(a, wd, f'apo_t_{c}', c).items():
                    free_t[(c,) + k] = v

            # FREE BINDER: delete every target chain -- the mirror image.
            b = gemmi.read_structure(cif); b.setup_entities()
            b.remove_ligands_and_waters(); b.setup_entities()
            for i in range(len(b[0]) - 1, -1, -1):
                if b[0][i].name != bnd: del b[0][i]
            b.setup_entities()
            free_b = propka_his_and_acids(b, wd, 'apo_b', bnd)

        # COMPOSE OVER EVERY TITRATABLE SITE ON BOTH PARTNERS -- HIS, ASP and GLU.
        #
        # The reviewer, 2026-10-05, named the defect: the old code read pKa values for HIS,
        # ASP and GLU alike, but only histidines were ever appended to `sites`, which left
        # every acid out of the multiplied product. They flagged that this is most damaging
        # precisely where a design's intervention is the introduction of an Asp or a Glu.
        #
        # They were right, and it was the worst possible place for that gap: the designed
        # intervention in this submission IS an acid introduction in most families (A22D,
        # S88D, S15D, L133E, T65D, S60D). propka already returned the acid pKa values in
        # both legs and this function discarded them, so the gate was blind to the very
        # residue each design was built around. Mechanism B -- a binder carboxylate reading
        # a target histidine -- cannot be seen by a histidine-only product at all.
        #
        # `product` is now the all-site product. `product_his_only` is kept alongside it as
        # the previous basis so the two can be compared directly rather than swapped
        # silently; that comparison is the pH sensitivity analysis, not a confidence interval.
        #
        # A site whose pKa is missing in either leg is now RECORDED in `unassessed` instead
        # of being skipped by a bare `continue`. A silent skip makes an unmeasured site
        # indistinguishable from an absent one, and makes the product's site count a lie.
        sites, unassessed = {}, {}

        def add(partner, label, rn, num, f, bo, chain):
            if f is None or bo is None:
                unassessed[f'{partner}:{label}'] = dict(
                    partner=partner, resname=rn, resnum=num,
                    free=None if f is None else round(f, 2),
                    bound=None if bo is None else round(bo, 2),
                    reason=('no pKa in either leg' if f is None and bo is None else
                            'no free-leg pKa' if f is None else 'no bound-leg pKa'))
                return
            sites[f'{partner}:{label}'] = dict(partner=partner, resname=rn, resnum=num,
                                               chain=chain,
                                               free=round(f, 2), bound=round(bo, 2),
                                               ratio=round(_link(f, bo), 4))

        # TARGET: histidines first, named from the construct's fingerprint, then its acids.
        # With more than one protomer the chain goes INTO the label, or A73 and B73 would
        # collide in the dict and one of them would silently vanish from the product.
        multi = len(tgts) > 1
        for c in tgts:
            pre = f'{c}/' if multi else ''
            for num, nm in names.items():
                add('target', f'{pre}{nm}', 'HIS', num,
                    free_t.get((c, 'HIS', num)), bound_t.get((c, 'HIS', num)), c)
            for (cc, rn, num) in sorted(k for k in set(free_t) | set(bound_t) if k[0] == c):
                if rn == 'HIS' and num in names: continue      # already added under its name
                add('target', f'{pre}{rn}{num}', rn, num,
                    free_t.get((c, rn, num)), bound_t.get((c, rn, num)), c)
        # BINDER: every titratable residue, not a curated list
        for (rn, num) in sorted(set(free_b) | set(bound_b)):
            add('binder', f'{rn}{num}', rn, num,
                free_b.get((rn, num)), bound_b.get((rn, num)), bnd)

        if not sites:
            return {'file': os.path.basename(cif), 'error': 'no titratable site measured',
                    'unassessed': unassessed, 'n_unassessed': len(unassessed)}

        # Guard 3 + 4, per site
        for k, v in sites.items():
            v['implied_pka_bound'] = implied_pka_bound(v['ratio'], v['free'])
            lo_p, hi_p = PKA_PLAUSIBLE_BY_RES.get(v['resname'], PKA_PLAUSIBLE)
            v['plausible_window'] = [lo_p, hi_p]
            v['implausible'] = (v['implied_pka_bound'] is None or
                                not (lo_p <= v['implied_pka_bound'] <= hi_p))
            # The other side may be SEVERAL chains: a binder histidine's nearest counter-charge
            # can sit on any protomer of a homotrimer, and taking only one would under-report
            # the partner and flag a real site `no_partner`.
            mine = v['chain']
            others = [bnd] if v['partner'] == 'target' else list(tgts)
            ds = [cross_partner_distance(st[0], mine, v['resnum'], v['resname'], o)
                  for o in others]
            ds = [x for x in ds if x is not None]
            d = min(ds) if ds else None
            v['partner_dist'] = d
            v['no_partner'] = (d is None or d > PARTNER_CUT)
            v['moved'] = abs(v['ratio'] - 1.0) >= NOISE

        # Guard 1: compose over EVERYTHING
        prod = 1.0
        for v in sites.values(): prod *= v['ratio']
        helpful = [v['ratio'] for v in sites.values() if v['ratio'] > 1.0]
        k_moved = sum(1 for v in sites.values() if v['moved'])
        # The previous basis, kept for a side-by-side comparison rather than replaced
        # silently. Any difference between these two IS the effect of the acids.
        his = [v['ratio'] for v in sites.values() if v['resname'] == 'HIS']
        prod_his = math.prod(his) if his else 1.0
        acids = [v for v in sites.values() if v['resname'] in ('ASP', 'GLU')]
        # Guard 2: ceiling
        ceiling = ONE_PROTON_BOUND ** max(k_moved, 1)
        verdict = 'IMPOSSIBLE' if prod > ceiling * 1.001 else 'ok'
        contributing = {k: v for k, v in sites.items() if v['moved']}
        return dict(file=os.path.basename(cif)[:-4], path=cif, target_chain='+'.join(tgts), binder_chain=bnd,
                    product=round(prod, 4),
                    product_his_only=round(prod_his, 4),
                    product_helpful_only_DO_NOT_USE=round(math.prod(helpful) if helpful else 1.0, 4),
                    n_sites=len(sites), n_his=len(his), n_acid=len(acids),
                    n_acid_moved=sum(1 for v in acids if v['moved']),
                    n_moved=k_moved, ceiling_for_n_moved=round(ceiling, 2),
                    verdict=verdict,
                    n_implausible=sum(1 for v in sites.values() if v['moved'] and v['implausible']),
                    n_no_partner=sum(1 for v in sites.values() if v['moved'] and v['no_partner']),
                    n_unassessed=len(unassessed), unassessed=unassessed,
                    contributing=contributing, sites=sites)
    except Exception as e:
        return {'file': os.path.basename(cif), 'error': f'{type(e).__name__}: {e}'}


def reconstruct(parent, singles, double):
    """Guard 5. Is the double's product the product of its single-mutant effects?

    effect(X) = product(X) / product(parent).  Expected double effect = prod(effects).
    Returns a dict including whether the observed effect is within tolerance.
    """
    if not parent or 'product' not in parent: return {'error': 'no parent product'}
    p0 = parent['product']
    if p0 <= 0: return {'error': 'parent product is zero'}
    eff = {}
    for nm, s in singles.items():
        if s and 'product' in s: eff[nm] = s['product'] / p0
    if not eff: return {'error': 'no single-mutant products'}
    expected = p0
    for e in eff.values(): expected *= e
    observed = double.get('product') if double else None
    if observed is None: return {'error': 'no double product'}
    ratio = observed / expected if expected else None
    return dict(parent_product=round(p0, 4),
                single_effects={k: round(v, 4) for k, v in eff.items()},
                expected_double=round(expected, 4), observed_double=round(observed, 4),
                observed_over_expected=round(ratio, 4) if ratio else None,
                additive_in_log=bool(ratio and 0.7 <= ratio <= 1.43),
                note=('reconstructs from its parts' if ratio and 0.7 <= ratio <= 1.43 else
                      'DOES NOT reconstruct -- epistasis or a bug, and this script '
                      'cannot distinguish them. Do not report the product as a result.'))


def selftest():
    # the physics
    # 40 / -20 stand in for the pKa_bound limits; 1e6 overflows the 10**x in link().
    assert abs(link(6.22, 40.0, 6.5, 7.4) - 5.554) < 0.01   # P1 reference
    assert abs(link(6.22, -20.0, 6.5, 7.4) - 0.699) < 0.01  # P1 reference
    assert abs(ONE_PROTON_BOUND - 7.943) < 0.001
    assert abs(link(6.0, 6.0, 6.5, 7.4) - 1.0) < 1e-9, 'no shift must give exactly 1.0'
    assert abs(link(6.0, 6.0, 6.0, 7.4) - 1.0) < 1e-9, 'true at either pH pair'
    # inversion round-trips
    for free, pb in [(6.22, 9.0), (6.22, 7.0), (4.93, 8.0), (6.5, 6.6)]:
        r = link(free, pb, 6.5, 7.4)
        assert abs(implied_pka_bound(r, free) - pb) < 0.05, (free, pb, r)
    assert implied_pka_bound(999.0, 6.22) is None, 'unreachable ratio must return None'
    # Guard 1: the historical bug. H370 6.056 x H383 1.503 while H433 reads 0.723
    sites = {'a': 6.056, 'b': 1.503, 'c': 0.723}
    allp = 1.0
    for v in sites.values(): allp *= v
    helpful = 6.056 * 1.503
    assert abs(helpful - 9.102) < 0.01, helpful
    assert abs(allp - 6.581) < 0.01, allp
    assert allp < helpful, 'all-sites must be below helpful-only on this case'
    assert allp < ONE_PROTON_BOUND, 'and below the single-proton ceiling'
    # Guard 2: a product above the k-site ceiling is impossible.
    # One site caps at 7.943x, two at 63.1x, three at 501x.
    assert abs(ONE_PROTON_BOUND ** 2 - 63.096) < 0.01, ONE_PROTON_BOUND ** 2
    assert 50.0 > ONE_PROTON_BOUND ** 1, 'a 50x product needs >1 site'
    assert 50.0 < ONE_PROTON_BOUND ** 2, 'and 2 sites suffice for it'
    assert 600.0 > ONE_PROTON_BOUND ** 2, 'a 600x product needs >2 sites'
    # Guard 5: reconstruction
    par = {'product': 1.0}
    rec = reconstruct(par, {'m1': {'product': 3.0}, 'm2': {'product': 2.0}}, {'product': 6.0})
    assert rec['additive_in_log'] and abs(rec['observed_over_expected'] - 1.0) < 1e-9, rec
    rec2 = reconstruct(par, {'m1': {'product': 3.0}, 'm2': {'product': 2.0}}, {'product': 20.0})
    assert not rec2['additive_in_log'], rec2
    # a double that is merely the better single does NOT reconstruct as multiplicative
    rec3 = reconstruct({'product': 0.70}, {'S88D': {'product': 4.57}, 'S60D': {'product': 1.11}},
                       {'product': 4.59})
    assert not rec3['additive_in_log'], rec3
    print('selftest OK')
    # --- 2026-10-08: the pH pair and the direction of merit ---
    # P1 and P2 run OPPOSITE ways. A 5x ratio is a success for problem 1 and a failure for
    # problem 2, so a direction-blind gate would hand problem 2 its own antithesis as a win.
    set_problem(1)
    assert (PH_LO, PH_HI) == (6.5, 7.4) and abs(ONE_PROTON_BOUND - 7.943) < 0.001
    assert merits(5.0) == 'helps' and merits(0.2) == 'wrong_direction', 'P1 wants acid-tightening'
    set_problem(2)
    assert (PH_LO, PH_HI) == (6.0, 7.4) and abs(ONE_PROTON_BOUND - 25.119) < 0.001
    assert merits(5.0) == 'wrong_direction' and merits(0.2) == 'helps', 'P2 wants acid-weakening'
    assert merits(1.0) == 'noise' and merits(1.0 + NOISE / 2) == 'noise'
    # MUTATION: the ceiling must move WITH the pair. If set_problem changed one and not the
    # other, a problem-2 product of 20x would be flagged IMPOSSIBLE against P1's 7.943.
    assert ONE_PROTON_BOUND > 20.0, 'P2 ceiling must admit 20x on one site'
    set_problem(1)
    assert ONE_PROTON_BOUND < 20.0, 'P1 ceiling must refuse 20x on one site'
    # MUTATION: _link must carry the module pair, not ph_gate_all's default
    set_problem(2)
    assert abs(_link(6.22, -20.0) - 0.401) < 0.01, _link(6.22, -20.0)
    set_problem(1)
    assert abs(_link(6.22, -20.0) - 0.699) < 0.01, _link(6.22, -20.0)
    # MUTATION: no bare link( may survive outside _link and these known-answer tests, or a
    # future call site could silently read the wrong pH pair -- the fault this closes.
    src = open(os.path.abspath(__file__)).read()
    body = src[:src.index('def selftest(')]
    bare = [l for l in body.splitlines()
            if re.search(r'(?<![_\w])link\(', l) and 'def _link' not in l
            and 'PH_LO, PH_HI)' not in l and not l.strip().startswith('#')]
    assert not bare, f'bare link( outside _link: {bare}'
    # --- orientation on N chains, 2026-10-08 ---
    # Stubs, because the only thing family() touches is residue names and numbers.
    class _R:
        def __init__(s, n, i): s.name, s.seqid = n, type('S', (), {'num': i})()
    class _C:
        def __init__(s, nm, hs): s.name, s._r = nm, [_R('HIS', i) for i in hs]
        def __iter__(s): return iter(s._r)
    from ph_gate_all import TNF, CROP
    tnf_hs, crop_hs = sorted(TNF), sorted(CROP)
    # problem 2: three protomers plus one binder
    tg, bd, nm, why = orient_multi([_C('A', tnf_hs), _C('B', tnf_hs),
                                    _C('C', tnf_hs), _C('D', [5, 40])])
    assert why is None and tg == ['A', 'B', 'C'] and bd == 'D' and nm is TNF, (tg, bd, why)
    # problem 1: one target plus one binder -- the case that must not change
    tg1, bd1, nm1, why1 = orient_multi([_C('A', [7]), _C('B', crop_hs)])
    assert why1 is None and tg1 == ['B'] and bd1 == 'A' and nm1 is CROP, (tg1, bd1, why1)
    # and with the chain order reversed, because this project's poses do both
    tg2, bd2, _, why2 = orient_multi([_C('A', crop_hs), _C('B', [7])])
    assert why2 is None and tg2 == ['A'] and bd2 == 'B', (tg2, bd2, why2)
    # MUTATION: no recognised family must REFUSE, not guess a target
    assert orient_multi([_C('A', [1, 2]), _C('B', [3, 4])])[3] == 'target family unrecognised'
    # MUTATION: two non-target chains must refuse -- a second binder, or a stray chain, would
    # otherwise be silently treated as part of the target side
    why3 = orient_multi([_C('A', tnf_hs), _C('B', [5]), _C('C', [6])])[3]
    assert 'need exactly one binder' in why3, why3
    # MUTATION: zero non-target chains (apo target alone) must refuse too
    assert 'need exactly one binder' in orient_multi([_C('A', tnf_hs)])[3]
    # MUTATION: the free-TARGET leg must delete the binder ONLY. Deleting every non-target
    # chain -- what the two-chain code did -- removes the sibling protomers, and H73 sits at
    # the inter-protomer seam, so its free pKa would come from a monomer that never exists.
    src = open(os.path.abspath(__file__)).read()
    assert "if a[0][i].name == bnd: del a[0][i]" in src, 'free-target leg deletes by binder name'
    assert "free-target leg lost a protomer" in src, 'the protomer count is not asserted'
    assert "if b[0][i].name != bnd: del b[0][i]" in src, 'free-binder leg must keep only binder'
    print(f'  orientation: 3 protomers + binder -> {tg}/{bd}; 1+1 still {tg1}/{bd1} either order')
    print(f'  MUTATION: unrecognised family, 2 binders and 0 binders all refuse')
    print(f'  MUTATION: free-target leg deletes the binder only, protomer count asserted')
    print(f'  pH pair + direction: P1 {PROBLEMS[1][0]}/{PROBLEMS[1][1]} wants ratio > 1, '
          f'P2 {PROBLEMS[2][0]}/{PROBLEMS[2][1]} wants ratio < 1')
    print(f'  one-proton ceiling moves with the pair: 7.943x (P1) / 25.119x (P2)')
    print(f'  MUTATION: no bare link( call sites remain outside _link')
    print(f'  single-site range at pKa_free 6.22 : {link(6.22,-20.0,6.5,7.4):.3f}x to {link(6.22,40.0,6.5,7.4):.3f}x  (problem 1)')
    print(f'  one-proton bound                   : {ONE_PROTON_BOUND:.3f}x')
    print(f'  historical bug reproduced          : helpful-only 9.102x vs all-sites 6.581x')
    print(f'  S88D+S60D does NOT reconstruct     : observed 4.59 vs expected '
          f"{rec3['expected_double']} (ratio {rec3['observed_over_expected']})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cifs', nargs='*')
    ap.add_argument('--dir'); ap.add_argument('--tsv'); ap.add_argument('--json')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--problem', type=int, default=1, choices=sorted(PROBLEMS),
                    help='1 = bind 6.5 / silent 7.4 (default, EGFR). '
                         '2 = bind 7.4 / silent 6.0 (TNF-alpha). Sets the pH pair, the '
                         'one-proton ceiling AND the direction of merit together.')
    a = ap.parse_args()
    if a.selftest: selftest(); return
    set_problem(a.problem)
    print(f"problem {PROBLEM}: pH {PH_LO} vs {PH_HI}, one-proton ceiling "
          f"{ONE_PROTON_BOUND:.3f}x per site")
    print(f"  {MERIT_TEXT}\n")
    cifs = a.cifs or (sorted(glob.glob(os.path.join(a.dir, '**', '*.cif'), recursive=True))
                      if a.dir else [])
    if not cifs: raise SystemExit(__doc__)
    rows = []
    for c in cifs:
        r = score_pose(c)
        rows.append(r)
        if 'error' in r:
            print(f"SKIP  {r['file']}: {r['error']}"); continue
        flag = (' ** ' + r['verdict']) if r['verdict'] != 'ok' else ''
        print(f"{r['file'][:52]:<54} product {r['product']:>8.3f}  "
              f"{r['n_moved']}/{r['n_sites']} sites moved  ceiling {r['ceiling_for_n_moved']:.1f}"
              f"{flag}")
        for k, v in sorted(r['contributing'].items(), key=lambda kv: -kv[1]['ratio']):
            w = []
            if v['implausible']: w.append(f"IMPLAUSIBLE pKa_bound={v['implied_pka_bound']}")
            if v['no_partner']: w.append(f"NO COUNTER-CHARGE within {PARTNER_CUT}A"
                                         + (f" (nearest {v['partner_dist']}A)" if v['partner_dist'] else ""))
            print(f"      {k:<22} {v['ratio']:>7.3f}x  {merits(v['ratio']):<16}"
                  f"free {v['free']:.2f} -> bound {v['bound']:.2f}  {'  '.join(w)}")
    if a.json: json.dump(rows, open(a.json, 'w'), indent=1); print(f"\nwrote {a.json}")
    if a.tsv:
        import csv as _csv
        with open(a.tsv, 'w', newline='') as fh:
            w = _csv.writer(fh, delimiter='\t')
            w.writerow(['file','product','n_sites','n_moved','ceiling','verdict',
                        'n_implausible','n_no_partner'])
            for r in rows:
                if 'error' in r: continue
                w.writerow([r['file'], r['product'], r['n_sites'], r['n_moved'],
                            r['ceiling_for_n_moved'], r['verdict'], r['n_implausible'],
                            r['n_no_partner']])
        print(f"wrote {a.tsv}")


if __name__ == '__main__':
    main()
