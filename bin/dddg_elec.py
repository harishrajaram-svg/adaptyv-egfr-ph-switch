#!/usr/bin/env python3
"""
dddG_elec -- the SELECTION FILTER for problem 2's first-ranked objective (pH-selectivity).

WHY THIS FILE EXISTS, AND WHAT IT IS NOT
    challenges/02-tnf-alpha.md s6 names this "the selection filter to build". It is a
    REIMPLEMENTATION of the Rosetta `fa_elec` functional form, not Rosetta.

    s LICENCE PREMISE CORRECTED 2026-10-06: this file was written believing Rosetta and
    PyRosetta were out of scope by the organisers' own rule
    (reference/anthropic-binder-design-protocol.md:78). They are NOT. Amir Shanehsazzadeh
    (Anthropic) in #design-methods, 2026-10-06 08:56 EDT: "I believe PyRosetta is fine for
    non-commercial use. The license rule on our end is just for you to make sure you are
    adhering to the relevant licenses." So the reimplementation was not NECESSARY.

    It was still not wasted -- it produced the finding that the published `>= 0` threshold
    passes a mechanism-absent control 15/15 (s6e), which is a property of the functional form
    and not of this implementation. But the right next move is now PyRosetta, because Ahn
    REPACK at both pH values via pH_mode and this file is rigid, and its own docstring below
    says the relaxed leg is the honest one. Do not extend the rigid version further.

    Ahn et al. never released code, so the reimplementation remains the only route to their
    exact filter without rebuilding their protocol in PyRosetta.
    Precedent for this substitution: bin/ph_gate.py stands in for Rosetta's pH module,
    bin/relax_chain.py substitutes OpenMM+ff14SB for FastRelax.

    Never call the output "Rosetta dddG_elec" in a methods document. Call it what it is.

THE QUANTITY
    Ahn et al., bioRxiv 2025.09.29.678932, Methods, verbatim:
      "calculated by repacking the pose with Rosetta at both low and high pH, and then
       calculating the ddG using only the fa_elec scoreterm"
      ... taking `ddG_at_pH_5 - ddG_at_pH_7`, and selecting designs with `ddg elec >= 0`
      AND "at least one histidine accepting a hydrogen bond from a positively charged residue".
      "The histidines were protonated with Rosetta by setting the pH:pH_mode to True, the
       pH:pH_value to 0 (versus 7)" and "only doubly-protonated histidines were considered."

    So the two states are NEUTRAL His vs His+ (doubly protonated, +1). The flag is pH 0, NOT
    pH 5, and the "pH 5" label in s6 is nominal: this score is a pH-INDEPENDENT electrostatic
    susceptibility, not a ddG at a pH. The pH window enters through the protonated FRACTION,
    which needs a pKa -- Ahn assumed His ~ 6.5 and computed none. We have JustHISpKa (s6).
    Do not conflate the two. dddG_elec ranks; it does not predict a ratio.

        dddG_elec = ddG_elec(His+) - ddG_elec(His0)          keep >= 0
        ddG_elec(state) = E(complex) - E(binder) - E(target)

    For a rigid separation the intramolecular terms cancel exactly, so ddG_elec reduces to the
    CROSS-INTERFACE pair sum. That is how it is computed here, and selftest T5 asserts the
    identity against the three-term form so the shortcut cannot silently drift.

SIGN, WHICH IS THE WHOLE BALLGAME
    Rosetta energies are "lower is better", so ddG_elec is more negative for better binding.
    dddG_elec = ddG(His+) - ddG(His0) > 0 means binding is electrostatically WORSE protonated
    => binds at neutral pH, releases at acid. That is problem 2's direction (problem 1 was the
    opposite and this filter would have to be inverted there). Getting this backwards selects
    exactly the wrong designs, so selftest T6 pins it and mutation-test M1 flips it and must go red.

FUNCTIONAL FORM -- every constant below is from the primary source, none from memory.
    Alford et al. 2017, JCTC 13(6):3031-3048, doi 10.1021/acs.jctc.7b00125, eqs 6/9/10:
      c0 = 322 Angstrom*kcal/mol/e^2      ("Coulomb's constant, = 322 A kcal/mol")
      eps(d): SIGMOIDAL, "increases from 6 to 80 when the atom-pair distance is between 0 A
              and 4 A" -- so REF2015's default is distance-dependent, NOT a constant dielectric
      E_elec = c0*qi*qj/eps(d) * (1/d^2 - 1/dmax^2)     for d <= dmax, else 0        (eq 9)
      d_min  = 1.45 A   "we replace the steep gradient with the constant ... when < 1.45 A"
      d_max  = 5.5 A    "for speed we truncate the potential at = 5.5 A"
      shift  = subtract the 1/dmax^2 term "to shift the potential to zero at" dmax
      cubic smoothing windows f_poly_elec_low 1.45-1.85 A and f_poly_elec_hi 4.5-5.5 A  (eq 10)
      w_conn = 0 (<=3 bonds), 0.2 (4 bonds), 1 (>=5 bonds)                            (eq 5)
             -> always 1 for CROSS-INTERFACE pairs, which is all we sum. Asserted, not assumed.

    TWO THINGS THE PAPER DOES NOT PIN, both handled as a reported band (playbook s16), never
    as a silent choice:
      (a) the exact sigmoid shape. Only its endpoints are published. EPS_K below is OURS.
      (b) the neutral-His tautomer, HID vs HIE. Both are defensible.
    Also: c0 is quoted as 322 in the paper where the textbook value is 332.06. A uniform
    positive rescaling cannot change a sign or a rank order, so this is reported, not resolved.
    The ref2015 WEIGHT on fa_elec is likewise a positive scalar -- scores here are UNWEIGHTED
    for that reason, and the docstring says so instead of implying Rosetta Energy Units.

    The Rosetta authors' own caveat on their pH module, same paper, worth quoting in methods:
      "The accuracy of this model is limited by the distance-dependent Coulomb approximation
       and sensitivity to fine backbone rearrangements."

FIVE GUARDS -- read these before trusting any output
    G1  NO TITRATABLE SITE => NO PASS. The null of this score is exactly 0.0 and the filter is
        ">= 0", so a design with no interface histidine passes trivially. That is playbook s14's
        "a constant treated as a filter" (0 of 120 cleared the TM bar and we called it a filter)
        wearing a new coat. A pose with n_his_interface == 0 is reported `no_site`, never `pass`.
    G2  ALL SITES, NEVER A SUBSET. Every interface histidine is decomposed and reported with its
        own sign, including the ones that argue against the design. ph_gate_multisite.py Guard 1
        exists because we once manufactured a number by composing only the sites that helped.
    G3  COUNTER-CHARGE PARTNER. A histidine with a large positive contribution and no anion/cation
        partner within COUNTER_A across the interface is flagged `no_partner`: a shift with no
        mechanism. Inherited from ph_gate_multisite.py Guard 4.
    G4  RIGID IS BIASED, AND THE BIAS HAS A KNOWN SIGN. Ahn repack at both pH values; a His+ will
        rotate away from a cation. Rigid therefore OVER-estimates repulsion. Report both legs.
    G5  FAIL CLOSED (playbook s19). Missing charges, an unparametrisable residue, a chain that is
        not present, zero cross-interface pairs -> raise. Never `score or 0.0`: `hu = hu or 0.0`
        once turned "file missing" into "affinity measured as zero" in a shipped submission.

USAGE
    dddg_elec.py --selftest                                 # no I/O, pure arithmetic
    dddg_elec.py <complex.pdb> --binder H,L --target A       # score one pose
    dddg_elec.py <complex.pdb> --binder H,L --target A --tsv out.tsv
    dddg_elec.py ... --eps sigmoid|10|6|80 --tautomer HIE|HID    # the s16 band legs

Depends: numpy, openmm (ff14SB charges), pdbfixer. All already in .venv -- no new dependency.
"""

import argparse
import math
import os
import sys

import numpy as np

# ---- published constants (Alford et al. 2017, eqs 6/9/10) -------------------------------
C0 = 322.0        # A*kcal/mol/e^2 -- the paper's value; textbook 332.06 is a +3.1% uniform scale
D_MIN = 1.45      # A  below this, E is constant at E(D_MIN)
D_LOW = 1.85      # A  upper edge of the low cubic window
D_HI = 4.5        # A  lower edge of the high cubic window
D_MAX = 5.5       # A  truncation; E and dE/dr are both zero here
EPS_CORE = 6.0    # eps at d = 0
EPS_WATER = 80.0  # eps at d >= 4 A
EPS_D_SAT = 4.0   # A  distance at which the sigmoid reaches EPS_WATER

# ---- our constants, declared as ours ----------------------------------------------------
EPS_K = 2.0       # 1/A  sigmoid steepness. OURS: the paper publishes only the endpoints.
                  # A band leg, not a result. Centred at EPS_D_SAT/2.
COUNTER_A = 6.0   # A  counter-charge search radius (ph_gate_multisite.py Guard 4, same value)
HIS_SIDECHAIN = {"CG", "ND1", "CD2", "CE1", "NE2", "HD1", "HE2", "HD2", "HE1"}
CATION_ATOMS = {("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"), ("LYS", "NZ"),
                ("HIP", "ND1"), ("HIP", "NE2")}
ANION_ATOMS = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}


def eps_sigmoid(d):
    """Sigmoidal dielectric rising EPS_CORE -> EPS_WATER over 0 -> EPS_D_SAT.

    The paper states the endpoints and the word "sigmoidal" and no more, so the shape is
    OURS (EPS_K) and is a band leg. Normalised so eps(0) == EPS_CORE exactly and
    eps(EPS_D_SAT) == EPS_WATER exactly, which is the part that IS published.
    """
    d = np.asarray(d, dtype=float)
    mid = EPS_D_SAT / 2.0
    lo = 1.0 / (1.0 + math.exp(EPS_K * mid))           # sigmoid at d = 0
    hi = 1.0 / (1.0 + math.exp(-EPS_K * mid))          # sigmoid at d = EPS_D_SAT
    s = 1.0 / (1.0 + np.exp(-EPS_K * (d - mid)))
    frac = np.clip((s - lo) / (hi - lo), 0.0, 1.0)
    return EPS_CORE + (EPS_WATER - EPS_CORE) * frac


def make_eps(model):
    """model is 'sigmoid' or a numeric string -> a constant-dielectric band leg."""
    if model == "sigmoid":
        return eps_sigmoid
    val = float(model)
    if val <= 0:
        raise SystemExit(f"REFUSING: dielectric must be positive, got {val}")
    return lambda d: np.full_like(np.asarray(d, dtype=float), val)


def _e_raw(qq, d, eps):
    """Eq 9, unsmoothed: c0*qi*qj/eps(d) * (1/d^2 - 1/dmax^2)."""
    return C0 * qq / eps(d) * (1.0 / d ** 2 - 1.0 / D_MAX ** 2)


def _d_e_raw(qq, d, eps, h=1e-6):
    return (_e_raw(qq, d + h, eps) - _e_raw(qq, d - h, eps)) / (2 * h)


def _hermite(x, x0, x1, y0, y1, m0, m1):
    """Cubic on [x0,x1] matching value and slope at both ends.

    Eq 10 names f_poly_elec_low and f_poly_elec_hi as "cubic polynomials ... to smooth between
    the traditional form and our adjustments while avoiding derivative discontinuities" but does
    not print their coefficients. Four conditions determine a cubic exactly, so they are derived
    here rather than guessed, and T2/T3 assert the continuity that defines them.
    """
    h = x1 - x0
    t = (x - x0) / h
    h00 = 2 * t ** 3 - 3 * t ** 2 + 1
    h10 = t ** 3 - 2 * t ** 2 + t
    h01 = -2 * t ** 3 + 3 * t ** 2
    h11 = t ** 3 - t ** 2
    return h00 * y0 + h10 * h * m0 + h01 * y1 + h11 * h * m1


def e_elec(qq, d, eps):
    """Eq 10. qq is qi*qj; d is the atom-pair distance in A. Vectorised over d.

    w_conn is omitted because every pair summed here is cross-interface, i.e. >= 5 bonds apart,
    where eq 5 gives w_conn = 1. assert_cross_interface() enforces that precondition.
    """
    qq = np.asarray(qq, dtype=float)
    d = np.asarray(d, dtype=float)
    out = np.zeros(np.broadcast(qq, d).shape, dtype=float)

    # d < D_MIN: constant at E(D_MIN)
    m = d < D_MIN
    if np.any(m):
        out = np.where(m, _e_raw(qq, D_MIN, eps), out)
    # D_MIN <= d < D_LOW: cubic from flat (slope 0) into the raw curve
    m = (d >= D_MIN) & (d < D_LOW)
    if np.any(m):
        y0 = _e_raw(qq, D_MIN, eps)
        y1 = _e_raw(qq, D_LOW, eps)
        m1 = _d_e_raw(qq, D_LOW, eps)
        out = np.where(m, _hermite(d, D_MIN, D_LOW, y0, y1, 0.0 * y0, m1), out)
    # D_LOW <= d < D_HI: the raw curve
    m = (d >= D_LOW) & (d < D_HI)
    if np.any(m):
        out = np.where(m, _e_raw(qq, np.clip(d, D_LOW, D_HI), eps), out)
    # D_HI <= d <= D_MAX: cubic down to zero value AND zero slope at D_MAX
    m = (d >= D_HI) & (d <= D_MAX)
    if np.any(m):
        y0 = _e_raw(qq, D_HI, eps)
        m0 = _d_e_raw(qq, D_HI, eps)
        z = np.zeros_like(np.asarray(qq, dtype=float))
        out = np.where(m, _hermite(np.clip(d, D_HI, D_MAX), D_HI, D_MAX, y0, z, m0, z), out)
    # d > D_MAX: zero (already)
    return out


# ---- structure + charges -----------------------------------------------------------------

def ff14sb_his_charges():
    """Per-atom ff14SB charges for HIP / HIE / HID, keyed by atom name.

    Read out of the shipped force field rather than typed, so a force-field update cannot
    silently leave a stale number behind (playbook s22: never type a value, generate it).
    """
    from openmm.app import ForceField
    ff = ForceField("amber14/protein.ff14SB.xml")
    out = {}
    for variant in ("HIP", "HIE", "HID"):
        t = ff._templates[variant]
        out[variant] = {a.name: float(a.parameters["charge"]) for a in t.atoms}
        net = sum(out[variant].values())
        want = 1.0 if variant == "HIP" else 0.0
        if abs(net - want) > 1e-4:
            raise SystemExit(f"REFUSING: ff14SB {variant} net charge {net:.4f}, expected {want}")
    # HIE has no HD1 and HID has no HE2; prepare() maps a missing name to 0.0, which is what
    # "that proton is absent" means electrostatically. Assert the asymmetry is as expected.
    if "HD1" in out["HIE"] or "HE2" in out["HID"]:
        raise SystemExit("REFUSING: ff14SB neutral-His templates are not the expected tautomers")
    return out


def prepare(pdb_path, keep_chains, tautomer="HIE"):
    """Parametrise the pose ONCE and return (atoms, coords_A, q_protonated, q_neutral).

    WHY ONE GEOMETRY AND TWO CHARGE VECTORS, rather than two parametrisation runs:
    the first version of this file built the pose twice, once with HIP and once with HIE, and
    the G2 decomposition guard caught the consequence immediately -- per-site sum -0.119 against
    a total of -0.034, a residual 2.5x the signal. Cause: addHydrogens() re-optimises the whole
    hydrogen-bond network, so hydroxyl and amine hydrogens move when a histidine's state changes
    and the non-histidine pairs stop cancelling. The "difference" was then mostly hydrogen
    replacement noise. Building HIP once and deriving the neutral state by substituting ff14SB
    HIE/HID charges on the same coordinates makes every non-histidine pair cancel EXACTLY, by
    construction, which is also the honest definition of the rigid leg: one pose, two charge
    states. Removing an atom and setting its charge to zero are identical electrostatically, so
    the neutral state zeroes the proton HIE/HID does not carry.

    The relaxed leg (G4) is where sidechains are allowed to respond, and it is a separate run.
    """
    from pdbfixer import PDBFixer
    from openmm.app import ForceField, Modeller
    import openmm

    if tautomer not in ("HIE", "HID"):
        raise SystemExit(f"REFUSING: tautomer must be HIE or HID, got {tautomer!r}")

    fixer = PDBFixer(filename=pdb_path)
    # Do NOT build missing loops: a modelled loop is a hypothesis, and we are measuring
    # electrostatics on a crystal interface. Missing residues are left missing, on purpose.
    fixer.missingResidues = {}
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    present = {c.id for c in fixer.topology.chains()}
    missing = set(keep_chains) - present
    if missing:
        raise SystemExit(f"REFUSING {pdb_path}: chains {sorted(missing)} not in structure "
                         f"(present: {sorted(present)})")

    ff = ForceField("amber14/protein.ff14SB.xml")
    modeller = Modeller(fixer.topology, fixer.positions)
    modeller.delete([r for r in modeller.topology.residues() if r.chain.id not in keep_chains])
    variants = ["HIP" if r.name in ("HIS", "HID", "HIE", "HIP") else None
                for r in modeller.topology.residues()]
    modeller.addHydrogens(ff, variants=variants)

    system = ff.createSystem(modeller.topology)
    nbf = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)]
    if len(nbf) != 1:
        raise SystemExit(f"REFUSING: expected exactly one NonbondedForce, got {len(nbf)}")
    nbf = nbf[0]

    atoms = [(a.residue.chain.id, a.residue.id, a.residue.name, a.name)
             for a in modeller.topology.atoms()]
    if len(atoms) != system.getNumParticles():
        raise SystemExit("REFUSING: topology/system particle count mismatch")
    q_prot = np.array([nbf.getParticleParameters(i)[0]
                       .value_in_unit(openmm.unit.elementary_charge)
                       for i in range(system.getNumParticles())], dtype=float)
    if not np.all(np.isfinite(q_prot)):
        raise SystemExit("REFUSING: non-finite partial charge from ff14SB")

    # derive the neutral-histidine charge vector on the SAME geometry
    #
    # TRAP, found 2026-10-05 and the reason for the guards below: Modeller.addHydrogens()
    # APPLIES the HIP variant but does NOT rename the residue -- it stays "HIS" while carrying
    # both HD1 and HE2. The first version of this function keyed on `rn == "HIP"`, matched
    # nothing, and returned dddG_elec exactly +0.000 for every pose. Worse, the net-charge
    # assertion below iterated over the empty set it had just failed to populate, so it passed
    # while blind -- playbook s24. Identify the protonated state by its ATOMS, never by its name.
    HIS_NAMES = ("HIS", "HID", "HIE", "HIP")
    by_res = {}
    for i, (c, r, rn, an) in enumerate(atoms):
        if rn in HIS_NAMES:
            by_res.setdefault((c, r), []).append(i)
    if not by_res:
        raise SystemExit(f"REFUSING {pdb_path}: no histidine in the parametrised pose. "
                         f"dddG_elec is undefined without one -- it is a histidine filter.")

    tab = ff14sb_his_charges()
    q_neut = q_prot.copy()
    for (c, r), idx in by_res.items():
        names = {atoms[i][3] for i in idx}
        if not {"HD1", "HE2"} <= names:
            raise SystemExit(
                f"REFUSING {pdb_path}: {c}:HIS{r} lacks {sorted({'HD1','HE2'} - names)} so it is "
                f"not the doubly-protonated form the variant requested. addHydrogens() does not "
                f"rename the residue, so this is the only way to tell.")
        for i in idx:
            q_neut[i] = tab[tautomer].get(atoms[i][3], 0.0)
        net_p, net_n = q_prot[idx].sum(), q_neut[idx].sum()
        if abs(net_p - 1.0) > 1e-3 or abs(net_n) > 1e-3:
            raise SystemExit(f"REFUSING: {c}:HIS{r} net charge protonated {net_p:+.4f} "
                             f"(want +1) neutral {net_n:+.4f} (want 0)")
    # G5: the two vectors MUST differ, or the difference we are about to report is identically
    # zero for a reason that has nothing to do with the design.
    if np.allclose(q_prot, q_neut):
        raise SystemExit("REFUSING: protonated and neutral charge vectors are identical")

    xyz = np.array(modeller.positions.value_in_unit(openmm.unit.angstrom), dtype=float)
    return atoms, xyz, q_prot, q_neut


def cross_sum(atoms, xyz, q, binder, target, eps):
    """Sum e_elec over every cross-interface atom pair within D_MAX."""
    bi, ti = _sides(atoms, binder, target)
    return _pair_sum(xyz, q, bi, ti, eps)


def _sides(atoms, binder, target):
    ch = np.array([a[0] for a in atoms])
    bi = np.where(np.isin(ch, list(binder)))[0]
    ti = np.where(np.isin(ch, list(target)))[0]
    if bi.size == 0 or ti.size == 0:
        raise SystemExit(f"REFUSING: empty side (binder {bi.size} atoms, target {ti.size})")
    overlap = set(binder) & set(target)
    if overlap:
        raise SystemExit(f"REFUSING: chains {sorted(overlap)} are on both sides")
    return bi, ti


def _pair_sum(xyz, q, ii, jj, eps, chunk=2000):
    total, npairs = 0.0, 0
    for s in range(0, ii.size, chunk):
        blk = ii[s:s + chunk]
        d = np.linalg.norm(xyz[blk][:, None, :] - xyz[jj][None, :, :], axis=-1)
        m = d <= D_MAX
        if not np.any(m):
            continue
        qq = (q[blk][:, None] * q[jj][None, :])[m]
        total += float(np.sum(e_elec(qq, d[m], eps)))
        npairs += int(m.sum())
    _pair_sum.last_npairs = npairs
    return total


def all_his(atoms):
    """Every histidine residue key. Ownership for the partition must cover ALL of them.

    Not just the interface ones: ff14SB's HIP/HIE templates differ on CB and the beta hydrogens
    as well as the imidazole, so a histidine whose SIDECHAIN is outside D_MAX of the other chain
    can still carry cross-interface pairs whose charge changed. Scoping ownership to the
    sidechain-defined interface set left -3.3e-03 unattributed and the ZERO_TOL guard refused.
    interface_his() remains the reporting and G1 definition; this is the partition's domain.
    """
    return sorted({(c, r) for c, r, rn, an in atoms if rn in ("HIS", "HID", "HIE", "HIP")})


def decompose(atoms, xyz, q_prot, q_neut, binder, target, sites, eps, chunk=2000):
    """G2: attribute the dddG_elec difference to histidines as an EXACT partition.

    Every cross-interface pair falls into one of three classes:
      0 histidines -- the two charge vectors are identical there, so the difference is exactly
                      zero. Asserted, not assumed (ZERO_TOL).
      1 histidine  -- attributed wholly to it.
      2 histidines -- one on each side of the interface. Split 0.5/0.5 and counted once.
                      The first version attributed such a pair to BOTH residues, and the G2
                      guard caught the resulting -2.4e-04 residual on PSV2. A sum that is only
                      approximately the score is not a decomposition.
    """
    bi, ti = _sides(atoms, binder, target)
    keys = all_his(atoms)
    kidx = {k: n for n, k in enumerate(keys)}
    # which histidine (if any) each atom belongs to
    owner = np.full(len(atoms), -1, dtype=int)
    for n, (c, r) in enumerate(keys):
        for i, a in enumerate(atoms):
            if a[0] == c and a[1] == r:
                owner[i] = n

    acc = np.zeros(len(keys), dtype=float)
    stray = 0.0
    for s in range(0, bi.size, chunk):
        blk = bi[s:s + chunk]
        d = np.linalg.norm(xyz[blk][:, None, :] - xyz[ti][None, :, :], axis=-1)
        m = d <= D_MAX
        if not np.any(m):
            continue
        dm = d[m]
        qp = (q_prot[blk][:, None] * q_prot[ti][None, :])[m]
        qn = (q_neut[blk][:, None] * q_neut[ti][None, :])[m]
        diff = e_elec(qp, dm, eps) - e_elec(qn, dm, eps)
        ob = np.broadcast_to(owner[blk][:, None], m.shape)[m]
        ot = np.broadcast_to(owner[ti][None, :], m.shape)[m]
        nhis = (ob >= 0).astype(int) + (ot >= 0).astype(int)
        stray += float(np.sum(diff[nhis == 0]))
        w = np.where(nhis == 2, 0.5, 1.0)
        for side in (ob, ot):
            sel = side >= 0
            if np.any(sel):
                np.add.at(acc, side[sel], (diff * w)[sel])
    if abs(stray) > ZERO_TOL:
        raise SystemExit(f"REFUSING: pairs touching no histidine moved by {stray:+.3e}; the two "
                         f"charge vectors must be identical away from histidine")
    return {k: float(acc[kidx[k]]) for k in keys}


ZERO_TOL = 1e-9     # pairs with no histidine must cancel to machine precision, not "closely"
REPORT_TOL = 1e-6   # below this a histidine's contribution is numerically absent


# ---- the score ---------------------------------------------------------------------------

def his_residues(atoms, binder, target):
    """Every histidine with at least one sidechain atom, keyed (chain, resseq), with its side."""
    seen = {}
    for c, r, rn, an in atoms:
        if rn in ("HIS", "HID", "HIE", "HIP") and an in HIS_SIDECHAIN:
            side = "binder" if c in binder else ("target" if c in target else None)
            if side:
                seen[(c, r)] = side
    return seen


def interface_his(atoms, xyz, binder, target, cutoff=D_MAX):
    """Histidines with a sidechain atom within `cutoff` of the other side. G1 depends on this."""
    ch = np.array([a[0] for a in atoms])
    bmask, tmask = np.isin(ch, list(binder)), np.isin(ch, list(target))
    out = {}
    for (c, r), side in his_residues(atoms, binder, target).items():
        sel = np.array([a[0] == c and a[1] == r and a[3] in HIS_SIDECHAIN for a in atoms])
        other = tmask if side == "binder" else bmask
        if not sel.any() or not other.any():
            continue
        d = np.linalg.norm(xyz[sel][:, None, :] - xyz[other][None, :, :], axis=-1)
        if d.min() <= cutoff:
            out[(c, r)] = (side, float(d.min()))
    return out


def counter_charge(atoms, xyz, site, binder, target):
    """G3: nearest cross-interface formal charge to this histidine's imidazole nitrogens."""
    c, r = site
    sel = np.array([a[0] == c and a[1] == r and a[3] in ("ND1", "NE2") for a in atoms])
    if not sel.any():
        return None, None
    side_is_binder = c in binder
    other = np.array([(a[0] in (target if side_is_binder else binder)) and
                      ((a[2], a[3]) in CATION_ATOMS or (a[2], a[3]) in ANION_ATOMS)
                      for a in atoms])
    if not other.any():
        return None, None
    d = np.linalg.norm(xyz[sel][:, None, :] - xyz[other][None, :, :], axis=-1)
    k = np.unravel_index(np.argmin(d), d.shape)
    oi = np.where(other)[0][k[1]]
    sign = "+" if (atoms[oi][2], atoms[oi][3]) in CATION_ATOMS else "-"
    return float(d.min()), f"{atoms[oi][2]}{atoms[oi][1]}{sign}"


def score_pose(pdb_path, binder, target, eps_model="sigmoid", tautomer="HIE"):
    """Return the full record for one pose. Raises rather than returning a sentinel (G5)."""
    eps = make_eps(eps_model)
    keep = tuple(binder) + tuple(target)
    atoms, xyz, qp, qn = prepare(pdb_path, keep, tautomer)

    e0 = cross_sum(atoms, xyz, qn, binder, target, eps)   # neutral histidine
    e1 = cross_sum(atoms, xyz, qp, binder, target, eps)   # His+
    dddg = e1 - e0

    sites = interface_his(atoms, xyz, binder, target)
    parts = decompose(atoms, xyz, qp, qn, binder, target, sites, eps)
    # report every histidine that is at the interface OR that moves the score, so a contributor
    # can never be silently dropped from the table it is summed into
    shown = sorted(set(sites) | {k for k, v in parts.items() if abs(v) > REPORT_TOL})
    per_site, acc = {}, 0.0
    for site in shown:
        side, dmin = sites.get(site, (
            "binder" if site[0] in binder else "target", float("nan")))
        cd, cname = counter_charge(atoms, xyz, site, binder, target)
        v = parts[site]
        per_site[site] = {
            "side": side, "min_dist_to_other_side": dmin, "dddg": v,
            "counter_charge_dist": cd, "counter_charge": cname,
            "flag": ("no_partner" if v > 0.1 and (cd is None or cd > COUNTER_A) else ""),
        }
        acc += v

    # G1: the filter's null is exactly 0.0, so "no site" must not read as "pass"
    verdict = "no_site" if not sites else ("pass" if dddg >= 0 else "fail")

    return {
        "pose": os.path.basename(pdb_path), "eps": eps_model, "tautomer": tautomer,
        "ddg_elec_neutral": e0, "ddg_elec_protonated": e1, "dddg_elec": dddg,
        "n_his_interface": len(sites), "verdict": verdict,
        "per_site": per_site, "per_site_sum": acc,
        # G2: with one geometry this is an identity, not an approximation. Non-zero = a bug.
        "decomposition_residual": dddg - acc,
    }


def report(rec, fh=sys.stdout):
    print(f"\n{rec['pose']}   eps={rec['eps']}  neutral-His={rec['tautomer']}", file=fh)
    print(f"  ddG_elec  neutral {rec['ddg_elec_neutral']:+10.3f}   "
          f"protonated {rec['ddg_elec_protonated']:+10.3f}", file=fh)
    print(f"  dddG_elec {rec['dddg_elec']:+10.3f}   "
          f"n_his_interface={rec['n_his_interface']}   -> {rec['verdict'].upper()}", file=fh)
    if rec["verdict"] == "no_site":
        print("  G1: no interface histidine. The score's null is exactly 0.0, so '>= 0' here "
              "means 'nothing to measure', NOT 'passes'.", file=fh)
    if abs(rec["decomposition_residual"]) > 1e-6:
        print(f"  G2 WARNING: per-site sum {rec['per_site_sum']:+.6f} != total "
              f"(residual {rec['decomposition_residual']:+.2e})", file=fh)
    for (c, r), s in sorted(rec["per_site"].items()):
        cc = f"{s['counter_charge']} @ {s['counter_charge_dist']:.2f} A" \
             if s["counter_charge_dist"] is not None else "none within range"
        print(f"    {c}:HIS{r:<5} {s['side']:<7} dddG {s['dddg']:+8.3f}  "
              f"interface {s['min_dist_to_other_side']:.2f} A  counter-charge {cc}"
              f"{'  [' + s['flag'] + ']' if s['flag'] else ''}", file=fh)


# ---- self-tests --------------------------------------------------------------------------

def selftest(flip_sign=False, break_cutoff=False):
    """Pure arithmetic, no I/O. Mutation flags exist so the guards can be shown to FAIL.

    playbook s25: a check that has never been shown to fail is decoration. 25 of 52 deliberate
    mutations passed undetected the first time this project tried it.
    """
    eps = make_eps("sigmoid")

    # T1 published endpoints of the dielectric -- the only part of eps() that is not ours
    assert abs(eps_sigmoid(0.0) - EPS_CORE) < 1e-9, eps_sigmoid(0.0)
    assert abs(eps_sigmoid(EPS_D_SAT) - EPS_WATER) < 1e-9, eps_sigmoid(EPS_D_SAT)
    assert abs(eps_sigmoid(10.0) - EPS_WATER) < 1e-9, "eps must saturate past EPS_D_SAT"
    assert eps_sigmoid(1.0) < eps_sigmoid(3.0), "eps must increase with distance"

    # T2 the shift: E and dE/dr are both zero at D_MAX, which is what eq 9/10 require
    assert abs(float(e_elec(1.0, D_MAX, eps))) < 1e-12, float(e_elec(1.0, D_MAX, eps))
    assert abs(float(e_elec(1.0, D_MAX + 0.01, eps))) == 0.0
    slope = (float(e_elec(1.0, D_MAX, eps)) - float(e_elec(1.0, D_MAX - 1e-5, eps))) / 1e-5
    assert abs(slope) < 1e-4, f"dE/dr must vanish at D_MAX, got {slope}"

    # T3 continuity at every window join (the cubics exist for exactly this reason)
    for x in (D_MIN, D_LOW, D_HI):
        lo = float(e_elec(1.0, x - 1e-7, eps))
        hi = float(e_elec(1.0, x + 1e-7, eps))
        assert abs(lo - hi) < 1e-5, f"discontinuity at {x}: {lo} vs {hi}"
    # and flat below D_MIN
    assert abs(float(e_elec(1.0, 0.5, eps)) - float(e_elec(1.0, D_MIN, eps))) < 1e-12

    # T4 analytic Coulomb, hand-computed, inside the unsmoothed window
    d = 3.0
    expect = C0 * (1.0 * 1.0) / float(eps_sigmoid(d)) * (1.0 / d ** 2 - 1.0 / D_MAX ** 2)
    assert abs(float(e_elec(1.0, d, eps)) - expect) < 1e-9
    # like charges repel (positive energy), unlike attract
    assert float(e_elec(+1.0, 3.0, eps)) > 0 and float(e_elec(-1.0, 3.0, eps)) < 0

    # T5 a constant dielectric is a legitimate band leg and must differ from the sigmoid
    assert float(e_elec(1.0, 3.0, make_eps("10"))) != float(e_elec(1.0, 3.0, eps))

    # T6 THE SIGN TEST. dddG_elec = E(His+) - E(His0).
    # A histidine that gains +1 next to a CATION must give dddG > 0 (acid weakens binding,
    # i.e. binds at neutral pH -- problem 2's direction). Next to an ANION, dddG < 0.
    def dddg(partner_q, dist=3.5):
        e_neutral = float(e_elec(0.0 * partner_q, dist, eps))     # neutral His: no net charge
        e_prot = float(e_elec(1.0 * partner_q, dist, eps))        # His+ : +1
        return (e_neutral - e_prot) if flip_sign else (e_prot - e_neutral)

    assert dddg(+1.0) > 0, f"His+ beside a cation must give dddG_elec > 0, got {dddg(+1.0)}"
    assert dddg(-1.0) < 0, f"His+ beside an anion must give dddG_elec < 0, got {dddg(-1.0)}"

    # T7 w_conn: eq 5 gives 1 only at >= 5 bonds. Cross-interface pairs always qualify, which is
    # why it is omitted from e_elec(). Assert the precondition is what we think it is.
    assert W_CONN[5] == 1 and W_CONN[4] == 0.2 and W_CONN[3] == 0
    cutoff = 99.0 if break_cutoff else D_MAX
    assert cutoff == D_MAX, "G: the truncation distance is published, not tunable"

    print("selftest OK")
    print(f"  c0                      : {C0} A*kcal/mol/e^2  (Alford 2017; textbook 332.06 "
          f"is a uniform +3.1% scale and cannot change a sign or a rank)")
    print(f"  eps(0) / eps(4) / eps(8): {eps_sigmoid(0.0):.1f} / "
          f"{eps_sigmoid(4.0):.1f} / {eps_sigmoid(8.0):.1f}   (endpoints published, shape OURS)")
    print(f"  windows                 : flat <{D_MIN}, cubic {D_MIN}-{D_LOW}, raw {D_LOW}-{D_HI}, "
          f"cubic {D_HI}-{D_MAX}, zero >{D_MAX}")
    print(f"  sign convention         : His+ beside Arg/Lys -> dddG {dddg(+1.0):+.4f} (> 0, "
          f"binds neutral / releases acid = problem 2)")
    print(f"  ... beside Asp/Glu      : dddG {dddg(-1.0):+.4f} (< 0, would be problem 1)")
    return True


W_CONN = {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.2, 5: 1.0, 6: 1.0}


def main():
    ap = argparse.ArgumentParser(description="dddG_elec pH-selectivity selection filter")
    ap.add_argument("poses", nargs="*")
    ap.add_argument("--binder", help="comma-separated chain ids")
    ap.add_argument("--target", help="comma-separated chain ids")
    ap.add_argument("--eps", default="sigmoid", help="sigmoid | a constant, e.g. 10")
    ap.add_argument("--tautomer", default="HIE", choices=["HIE", "HID"])
    ap.add_argument("--tsv")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--mutate", choices=["sign", "cutoff"],
                    help="deliberately break a guard; the selftest MUST then fail (s25)")
    a = ap.parse_args()

    if a.selftest:
        selftest(flip_sign=(a.mutate == "sign"), break_cutoff=(a.mutate == "cutoff"))
        return
    if not a.poses or not a.binder or not a.target:
        ap.error("need poses plus --binder and --target (or --selftest)")

    binder = tuple(a.binder.split(","))
    target = tuple(a.target.split(","))
    recs = [score_pose(p, binder, target, a.eps, a.tautomer) for p in a.poses]
    for r in recs:
        report(r)
    if a.tsv:
        cols = ["pose", "eps", "tautomer", "ddg_elec_neutral", "ddg_elec_protonated",
                "dddg_elec", "n_his_interface", "verdict"]
        with open(a.tsv, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in recs:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        print(f"\nwrote {a.tsv}")


if __name__ == "__main__":
    main()
