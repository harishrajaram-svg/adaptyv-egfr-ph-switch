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
        "iptm_repred": 0.39, "plddt_binder_repred": 0.78,
        "frac_V": 0.05, "frac_G": 0.04, "n_C": 0,
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
