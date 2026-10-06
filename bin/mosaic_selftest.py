#!/usr/bin/env python3
"""Self-test for modal_mosaic.py step 4 -- no GPU, no Modal account, no money.

Everything here is the part of step 4 that a GPU cannot check for you: the target
preparation, the numbering arithmetic, and the geometry measurement. Three numbering
schemes are live in this project (ECD-positional, domain-III-local, BoltzGen-output) and
mixing them has already cost real money, so the indices get asserted against an
INDEPENDENT fact -- the TKQHGQF motif -- rather than against arithmetic that could be
wrong in the same direction twice.

usage: python3 bin/mosaic_selftest.py        (run from the repo root)
"""
import sys
import types
from pathlib import Path

TGT = Path("targets/egfr/egfr_ecd_6aru_renum.pdb")
MOD = Path("biomodals/modal_mosaic.py")


def load_module():
    """Import modal_mosaic.py with `modal` stubbed out, so no account is needed."""
    m = types.ModuleType("modal")

    class _Deco:
        def __call__(self, fn):
            return fn

    class _App:
        def __init__(self, *a, **k):
            pass

        def function(self, *a, **k):
            return _Deco()

        def local_entrypoint(self, *a, **k):
            return _Deco()

    class _Volume:
        @staticmethod
        def from_name(*a, **k):
            return _Volume()

        def commit(self):
            pass

    class _Image:
        def __getattr__(self, _):
            return lambda *a, **k: self

        @staticmethod
        def debian_slim(*a, **k):
            return _Image()

    m.App, m.Volume, m.Image = _App, _Volume, _Image
    sys.modules["modal"] = m
    ns = {"__name__": "modal_mosaic_selftest"}
    exec(compile(MOD.read_text(), str(MOD), "exec"), ns)
    return ns


def main():
    if not TGT.is_file():
        sys.exit(f"run me from the repo root: {TGT} not found")
    ns = load_module()
    prepare, geometry = ns["_prepare_target"], ns["_geometry"]

    # 1. the full ECD: the scheme everything else is expressed in
    _, seq, nums = prepare(TGT, "A", None)
    assert len(seq) == 609, len(seq)
    assert nums == list(range(1, 610)), "ECD renum is not contiguous 1..609"
    assert seq[405] == "H", f"positional 406 is {seq[405]}, not H"
    assert seq[402:409] == "TKQHGQF", seq[402:409]
    print("full ECD      609 res, H433 at positional 406, motif TKQHGQF")

    # 2. the default crop, and the exact indices the GPU run will use
    chain, seq, nums = prepare(TGT, "A", (284, 453))
    idx = {n: i for i, n in enumerate(nums)}
    epitope = [idx[p] for p in range(403, 410)]
    assert (len(seq), len(chain), nums[0], nums[-1]) == (170, 170, 284, 453)
    assert "".join(seq[i] for i in epitope) == "TKQHGQF", "epitope indices are wrong"
    assert idx[406] == 122, idx[406]
    # cross-check against the OTHER live scheme: targets/egfr/egfr_d3_6aru_renum.pdb is
    # canonical 311..480 renumbered 1..170, so H433 sits at 123 there, i.e. 0-based 122.
    print(f"crop 284:453  170 res, epitope idx {epitope}, anchor idx {idx[406]}"
          " (== domain-III-local 123, independent)")

    # 3. the geometry check, on planted carboxylates at known distances
    import gemmi

    st = gemmi.read_structure(str(TGT))
    st.setup_entities()
    st.remove_ligands_and_waters()
    st.remove_alternative_conformations()
    his = next(r for r in st[0]["A"] if r.seqid.num == 406)
    nd1 = next(a for a in his if a.name == "ND1")

    def complex_with(offsets):
        """Binder chain A of Asp residues whose carboxylate O sits `offset` A from ND1."""
        fake = gemmi.Structure()
        model = gemmi.Model("1")
        binder = gemmi.Chain("A")
        for n, (offset, oname) in enumerate(offsets, start=1):
            res = gemmi.Residue()
            res.name, res.seqid = "ASP", gemmi.SeqId(n, " ")
            for name, d in ((oname, offset), ("CA", offset + 2.0)):
                atom = gemmi.Atom()
                atom.name, atom.element = name, gemmi.Element(name[0])
                atom.pos = gemmi.Position(nd1.pos.x + d, nd1.pos.y, nd1.pos.z)
                res.add_atom(atom)
            binder.add_residue(res)
        target = gemmi.Chain("B")
        target.add_residue(his)
        model.add_chain(binder)
        model.add_chain(target)
        fake.add_model(model)
        return fake

    g = geometry(complex_with([(3.0, "OD1"), (9.0, "OD2")]), 406)
    assert g["acid_O_to_his_N"] == 3.0, g          # nearest O wins, not the nearest residue
    assert g["closest_acid"] == "ASP1", g
    assert g["geometry_pass"] is True, g
    assert geometry(complex_with([(3.0, "OD1")]), 999) == {}, \
        "a missing anchor must report nothing, not a number"
    print(f"geometry      {g}")

    # 3b. PROBLEM 2, MULTI-CHAIN. _geometry_p2 must find an anchor on ANY target protomer, not
    # only the first. The original _geometry read st[0][1] alone; on a trimer that silently
    # measured against protomer A whatever chain the anchor was declared on -- a wrong number,
    # not an error. This is also the only multi-chain coverage in this file: the single-chain
    # path passed clean while the trimer run died on a NameError 20 GPU-minutes in, because
    # nothing here exercised it.
    geometry_p2 = ns["_geometry_p2"]

    def trimer_with(his_offset, anchor_on_chain):
        """Binder + THREE target protomers; the Arg anchor sits on `anchor_on_chain` (0,1,2)."""
        fake = gemmi.Structure()
        model = gemmi.Model("1")
        binder = gemmi.Chain("A")
        res = gemmi.Residue()
        res.name, res.seqid = "HIS", gemmi.SeqId(1, " ")
        for name, d in (("ND1", his_offset), ("CA", his_offset + 2.0)):
            atom = gemmi.Atom()
            atom.name, atom.element = name, gemmi.Element("N" if name[0] == "N" else "C")
            atom.pos = gemmi.Position(nd1.pos.x + d, nd1.pos.y, nd1.pos.z)
            res.add_atom(atom)
        binder.add_residue(res)
        model.add_chain(binder)
        # EVERY protomer is ARG 27, because TNF is a HOMOTRIMER and that is what the real
        # renumbered target looks like. The previous version named the decoys ALA, which let
        # the ("ARG","LYS") filter disambiguate the search -- so the test asserted a property
        # the production input does not have, and passed while _geometry_p2 was measuring the
        # wrong protomer on every real run. A test whose fixture is easier than production is
        # worse than no test, because it is reported as coverage.
        for ci in range(3):
            tgt = gemmi.Chain("BCD"[ci])
            arg = gemmi.Residue()
            arg.name = "ARG"
            arg.seqid = gemmi.SeqId(27, " ")
            atom = gemmi.Atom()
            atom.name, atom.element = "NH2", gemmi.Element("N")
            # the decoys sit 60 A away, so measuring the wrong one is a visible failure
            shift = 0.0 if ci == anchor_on_chain else 60.0
            atom.pos = gemmi.Position(nd1.pos.x + shift, nd1.pos.y, nd1.pos.z)
            arg.add_atom(atom)
            tgt.add_residue(arg)
            model.add_chain(tgt)
        fake.add_model(model)
        return fake

    # The anchor index must SELECT the protomer, not search for it. Passing the index of the
    # protomer the loss was aimed at must give 3.0 A; passing any other index must give ~60 A.
    for which in (0, 1, 2):
        g2 = geometry_p2(trimer_with(3.0, which), 27, which)
        assert g2.get("his_N_to_cation_N") == 3.0, (which, g2)
        assert g2["geometry_pass"] is True, (which, g2)
        # MUTATION TEST: the wrong index must NOT quietly return the right answer. This is the
        # assertion the old fixture could not make, and the bug it would have caught.
        for other in (0, 1, 2):
            if other == which:
                continue
            gw = geometry_p2(trimer_with(3.0, which), 27, other)
            assert gw["geometry_pass"] is False, (
                f"anchor on protomer {which} measured from index {other} reported "
                f"{gw.get('his_N_to_cation_N')} -- the index is being ignored")
    far = geometry_p2(trimer_with(7.0, 2), 27, 2)
    assert far["geometry_pass"] is False, far
    assert geometry_p2(trimer_with(3.0, 1), 999, 1) == {}, \
        "a missing anchor must report nothing, not a number"
    assert geometry_p2(trimer_with(3.0, 1), 27, 9) == {}, \
        "an out-of-range anchor index must report nothing, not fall back to protomer 0"
    print("p2 trimer     all three protomers are ARG 27 (as in production); the anchor INDEX "
          "selects,\n              and every wrong index reports a miss; 7.0 A -> False; "
          "missing anchor and\n              out-of-range index both -> {}")

    # 3c. THE REDUCTION. HisNearCation first wrote `score.sum()`, and summing rewards total
    # histidine MASS near the cation rather than one histidine PLACED: the 5-step trimer smoke
    # came back with 12 histidines in 76 residues (15.8% vs ~2.3% natural) while `his_best` sat
    # at exactly 0.00 every step. This is a MUTATION TEST -- it asserts that the discarded `sum`
    # formulation FAILS, so the test cannot pass if someone reverts the reduction.
    import math

    # No numpy. bin/design-mosaic.sh runs this file with bare `python3`, and gate_sweep.py:32
    # records what hardcoding .venv/bin/python cost: two gates died with FileNotFoundError in
    # every fresh clone. his_reduce needs sort, slice and sum and nothing else, which is the
    # whole reason it takes its array module as an argument -- so a six-line backend over
    # lists exercises the REAL reduction and keeps this file runnable on any interpreter.
    class _Arr(list):
        def sum(self):
            return sum(self)

        def __getitem__(self, k):
            got = list.__getitem__(self, k)
            return _Arr(got) if isinstance(k, slice) else got

    np = type("np", (), {"array": staticmethod(_Arr),
                         "sort": staticmethod(lambda x: _Arr(sorted(x)))})

    his_reduce = ns["his_reduce"]
    d0, width, eps = 6.5, 1.5, 1e-3

    def loss(n_his, dist, reduce):
        """n_his histidines, all `dist` A from the cation, padded out to 76 positions."""
        sig = 1.0 / (1.0 + math.exp(-(d0 - dist) / width))
        score = np.array([sig] * n_his + [0.0] * (76 - n_his))
        return -math.log(float(reduce(score)) + eps)

    placed = (1, 3.0)        # one histidine where we want it
    spam = (12, 9.0)         # twelve loitering at a distance that is useless on its own
    top2 = lambda s: his_reduce(s, 2, np)
    bad_sum = lambda s: s.sum()

    assert loss(*spam, bad_sum) < loss(*placed, bad_sum), \
        "the sum formulation is supposed to prefer the spam -- if it does not, this test is " \
        "no longer pinning down the bug it was written for"
    assert loss(*placed, top2) < loss(*spam, top2), \
        f"top-2 must prefer the placed His: {loss(*placed, top2)} vs {loss(*spam, top2)}"
    # and k=2 must actually VALUE the second site -- s24 measured 0.40 pKa units per site, so
    # one site cannot reach the spec and a max() reduction would pay for only one.
    assert loss(2, 3.0, top2) < loss(1, 3.0, top2), "k=2 does not reward a second placed His"
    assert loss(3, 3.0, top2) == loss(2, 3.0, top2), "a THIRD His must buy nothing"
    print(f"reduction      sum prefers spam ({loss(*spam, bad_sum):+.3f} < "
          f"{loss(*placed, bad_sum):+.3f}) -- top2 prefers placed "
          f"({loss(*placed, top2):+.3f} < {loss(*spam, top2):+.3f}); 2nd His pays, 3rd is free")

    # 3d. THE SPECIES LEG's index mapping, on the REAL structures. Mouse epitope positions
    # are NOT the human ones -- a deletion near human positional 68 shifts the register by -3
    # at the anchor and -4 further along -- so this is checked against residue IDENTITY, not
    # against arithmetic that could be wrong in the same direction twice. If the mapping slips,
    # the mouse leg aims the mechanism at whatever sits at the stale index and still reports a
    # plausible number.
    prep_set = ns["_prepare_target_set"]
    spec_file = Path("analysis/02-tnf/epitope_conserved.json")
    mouse_pdb = Path("targets/tnf/tnf_mouse_trimer_renum.pdb")
    if spec_file.is_file() and mouse_pdb.is_file():
        import json

        spec = json.loads(spec_file.read_text())
        hu_epi, mo_epi = spec["human_positional"], spec["mouse_positional"]
        assert len(hu_epi) == len(mo_epi), (hu_epi, mo_epi)

        hu = prep_set("t-human", Path("targets/tnf/tnf_trimer_renum.pdb"), "A,B,C", None,
                      hu_epi, spec["anchor"]["human_positional"], "B", "his_near_cation")
        # A,B -- a DIMER, forced: a mouse trimer is 976 tokens and OOMs on an L40S
        # (34.24 GiB allocation, 48 GB card). See the note on target2_chain.
        mo = prep_set("t-mouse", mouse_pdb, "A,B", None, mo_epi,
                      spec["anchor"]["mouse_positional"], "A", "his_near_cation")

        # the anchor must be the SAME residue in both, or one pH term cannot serve both legs
        assert hu["anchor_aa"] == mo["anchor_aa"] == "R", (hu["anchor_aa"], mo["anchor_aa"])

        # and every mapped epitope position must carry the SAME residue in both species --
        # that is what "conserved" has to mean for the two legs to be the same objective
        def residues(prep, positions):
            _cid, _ch, seq, idx, _rn = prep["chains_prepared"][0]
            return "".join(seq[idx[p]] for p in positions)

        hu_res, mo_res = residues(hu, hu_epi), residues(mo, mo_epi)
        assert hu_res == mo_res, f"epitope residues differ: human {hu_res} vs mouse {mo_res}"

        # index counts must scale with protomer count, not silently collapse to one
        assert len(hu["epitope_idx"]) == len(hu_epi) * 3, len(hu["epitope_idx"])
        assert len(mo["epitope_idx"]) == len(mo_epi) * 2, len(mo["epitope_idx"])

        # MUTATION TEST: the human numbers must NOT work on the mouse target. If they did,
        # the mapping would be decoration and a stale copy would go unnoticed.
        try:
            bad = prep_set("t-bad", mouse_pdb, "A,B", None, hu_epi,
                           spec["anchor"]["human_positional"], "A", "his_near_cation")
            assert residues(bad, hu_epi) != hu_res, (
                "human epitope numbers reproduce the human residues on the MOUSE target -- "
                "the -3/-4 register shift is not being exercised, so this test is vacuous")
        except SystemExit:
            pass        # refusing outright is also a correct answer

        print(f"species leg   anchor R in both; {len(hu_epi)} mapped positions carry "
              f"{hu_res} in human and mouse alike; {len(hu['epitope_idx'])} vs "
              f"{len(mo['epitope_idx'])} target idx; human numbers do NOT work on mouse")
    else:
        raise SystemExit("species-leg fixtures missing: "
                         f"{spec_file} and {mouse_pdb} are required")

    # 4. a distant carboxylate must FAIL the 4.0 A bar -- the bar is the whole point.
    # The reported distance is the min over ND1 AND NE2, so it is not the planted offset;
    # the bar is what is being tested here, not the arithmetic.
    g = geometry(complex_with([(7.0, "OD1")]), 406)
    assert g["acid_O_to_his_N"] > 4.0 and g["geometry_pass"] is False, g
    print(f"geometry bar  {g['acid_O_to_his_N']} A -> geometry_pass="
          f"{g['geometry_pass']}, as it must")

    # 5. the metric audit, against the two bugs this file actually shipped on 2026-10-03
    check_row, loss_term_names, audit_c1 = (
        ns["check_row"], ns["loss_term_names"], ns["audit_c1"]
    )
    clean = {
        "design": "t", "length": 76, "seed": 0, "sequence": "A" * 76,
        "iptm_design": 0.42, "plddt_binder_design": 0.81,
        "iptm_repred": 0.39, "iptm_repred_best_pair": 0.41,
        "plddt_binder_repred": 0.78,
        "frac_V": 0.05, "frac_G": 0.04, "frac_H": 0.03, "n_C": 0,
        "acid_O_to_his_N": 3.1, "acid_CA_to_his_N": 5.4,
        "closest_acid": "ASP34", "geometry_pass": True,
    }
    assert check_row(clean, 76) == [], check_row(clean, 76)

    def with_(**kw):
        r = dict(clean)
        r.update(kw)
        return r

    cases = [
        # the real bug: IPTMLoss handed a zero-length sequence returns exactly 0.0
        ("iptm_repred exactly 0.0", with_(iptm_repred=0.0), "iptm_repred"),
        # the real bug: a distance reported as negative
        ("negative distance", with_(acid_CA_to_his_N=-39.59), "acid_CA_to_his_N"),
        ("NaN metric", with_(iptm_design=float("nan")), "NaN"),
        ("probability over 1", with_(plddt_binder_repred=1.4), "plddt_binder_repred"),
        ("frac_H impossible", with_(frac_H=1.3), "frac_H"),
        # the problem-2 distance columns, which METRIC_BOUNDS used to omit entirely
        ("p2 distance negative", with_(his_N_to_cation_N=-3.0), "his_N_to_cation_N"),
        ("p2 distance NaN", with_(his_CA_to_cation_N=float("nan")), "NaN"),
        ("cysteine present", with_(n_C=5), "n_C"),
        ("wrong sequence length", with_(sequence="A" * 70), "sequence is 70"),
        ("pass disagrees with distance",
         with_(acid_O_to_his_N=7.2, geometry_pass=True), "disagrees"),
        ("distance without a residue",
         with_(closest_acid=None), "present or absent together"),
    ]
    for label, row, needle in cases:
        found = check_row(row, 76)
        assert found, f"{label}: audit found nothing"
        assert any(needle in f for f in found), f"{label}: {found}"
    # a BAD DESIGN must still pass -- the audit judges measurement, not quality
    bad_design = with_(iptm_repred=0.01, acid_O_to_his_N=19.4, geometry_pass=False,
                       frac_G=0.33)
    assert check_row(bad_design, 76) == [], check_row(bad_design, 76)
    print(f"metric audit   clean row OK, {len(cases)} impossible rows all caught, "
          "a merely-bad design still passes")

    # 5b. D-P2-1's iptm gate, and the guard that was BLIND on problem 2.
    #
    # Two things are tested. First the gate: below iptm_repred 0.45 the verdict is withdrawn to
    # None ("n/a"), never False, because False asserts the residue was measured and found
    # misplaced when the binder was not bound at all. Second, a MUTATION TEST that the
    # pass-vs-distance guard now fires on his_near_cation rows -- it read only
    # `acid_O_to_his_N`, problem 1's key, so on every problem 2 row it found None and skipped,
    # passing while completely blind (§25).
    gate = ns["gate_geometry"]
    BAR = ns["GEOMETRY_IPTM_BAR"]
    assert BAR == 0.45, BAR

    def p2row(**kw):
        r = {"design": "t", "length": 76, "seed": 0, "sequence": "A" * 76,
             "iptm_design": 0.42, "plddt_binder_design": 0.81,
             "iptm_repred": 0.80, "plddt_binder_repred": 0.78,
             "frac_V": 0.05, "frac_G": 0.04, "frac_H": 0.03, "n_C": 0,
             "his_N_to_cation_N": 3.1, "his_CA_to_cation_N": 5.4,
             "closest_his": "HIS21", "geometry_pass": True}
        r.update(kw)
        return r

    # above the bar: verdict stands, and a clean p2 row passes the audit
    hi = gate(p2row())
    assert hi["geometry_pass"] is True and hi["geometry_read"] is True, hi
    assert check_row(hi, 76) == [], check_row(hi, 76)

    # below the bar: verdict withdrawn, DISTANCE KEPT, reason recorded
    lo = gate(p2row(iptm_repred=0.1776, his_N_to_cation_N=20.36, geometry_pass=False))
    assert lo["geometry_pass"] is None, lo
    assert lo["geometry_read"] is False, lo
    assert lo["his_N_to_cation_N"] == 20.36, "the measured distance must survive the gate"
    assert "0.45" in lo["geometry_gate"], lo["geometry_gate"]
    assert check_row(lo, 76) == [], check_row(lo, 76)

    # a 12-histidine 3.15 A "pass" below the bar is withdrawn too -- that is p2trimer02, the
    # run whose pass came from histidine density rather than placement
    sp = gate(p2row(iptm_repred=0.3123, his_N_to_cation_N=3.15, geometry_pass=True))
    assert sp["geometry_pass"] is None, "a LUCKY pass below the bar must also be withdrawn"

    # missing / NaN iptm must fail closed, not sail through
    for bad_iptm in (None, float("nan")):
        g = gate(p2row(iptm_repred=bad_iptm))
        assert g["geometry_pass"] is None and g["geometry_read"] is False, (bad_iptm, g)

    # #8: geometry simply NOT COMPUTED above the bar must not be flagged. gate_geometry
    # leaves geometry_pass ABSENT when _geometry* returns {}, and check_row used row.get(),
    # collapsing absent and present-and-None -- so a design with fine metrics was marked
    # METRICS SUSPECT and dropped from ranking by its own audit.
    nogeom = {k: v for k, v in p2row(iptm_repred=0.80).items()
              if k not in ("his_N_to_cation_N", "his_CA_to_cation_N",
                           "closest_his", "geometry_pass")}
    ng = gate(nogeom)
    assert ng["geometry_read"] is True, ng
    assert check_row(ng, 76) == [], (
        "a row whose geometry was never computed, above the bar, must not be flagged: "
        + repr(check_row(ng, 76)))

    # #10: topk=0 silently restores the summed form this function exists to replace
    try:
        his_reduce(np.array([0.1, 0.2, 0.3]), 0, np)
        raise AssertionError("his_reduce(topk=0) must raise, not sum the whole array")
    except ValueError:
        pass

    # MUTATION TEST: the p2 disagreement guard must now FIRE. Above the bar so the gate does
    # not withdraw the verdict, with a verdict that contradicts its own distance.
    blind = gate(p2row(iptm_repred=0.80, his_N_to_cation_N=19.9, geometry_pass=True))
    found = check_row(blind, 76)
    assert any("his_N_to_cation_N" in f and "disagrees" in f for f in found), (
        "the pass-vs-distance guard is still blind on his_near_cation rows: " + repr(found))
    # and the paired-presence rule too
    orphan = gate(p2row(closest_his=None))
    assert any("present or absent together" in f for f in check_row(orphan, 76)), \
        check_row(orphan, 76)
    # problem 1 rows must be unaffected by all of this
    p1 = gate({**clean, "iptm_repred": 0.80})
    assert p1["geometry_pass"] is True and check_row(p1, 76) == [], check_row(p1, 76)

    print(f"D-P2-1 gate    bar {BAR}; below it the verdict is n/a not False and the distance "
          f"survives;\n               a lucky 3.15 A pass is withdrawn; missing/NaN iptm fails "
          f"closed;\n               p2 disagreement guard now FIRES (was blind); p1 unaffected")

    # 6. C1 derived from a loss tree, not asserted
    class Combo:
        def __init__(self, *members):
            self.l = list(members)

    def term(name, inner=None):
        t = type(name, (), {})()
        t.loss = inner
        return t

    tree = Combo(
        term("BinderTargetIPTM"), term("BinderTargetPAE"),
        term("BinderTargetContact"), term("PLDDTLoss"),
        term("ESMFoldGlobularity"),           # model-agnostic, must be ALLOWED
        term("AcidNearHis"),
    )
    wrapped = term("NoCys", inner=term("Boltz2Loss", inner=tree))
    names = loss_term_names(wrapped)
    assert "AcidNearHis" in names and "ESMFoldGlobularity" in names, names
    assert audit_c1(names, ["Boltz2"]) == [], audit_c1(names, ["Boltz2"])
    # and it must FAIL on the two things C1 forbids
    bad_terms = loss_term_names(Combo(term("IPSAE_min"), term("PLDDTLoss")))
    assert any("ipSAE" in v for v in audit_c1(bad_terms, ["Boltz2"])), bad_terms
    assert any("ESMFold2" in v for v in audit_c1(names, ["ESMFold2"]))
    print(f"C1 audit       {len(names)} terms walked, ESMFoldGlobularity allowed, "
          "IPSAE_min and an ESMFold2 model both rejected")

    print("\nselftest OK")


if __name__ == "__main__":
    main()
