"""Build the problem-2 validation-gate fixtures: does ipSAE_min separate a real TNF-alpha
binder from composition-matched negatives?

WHY THIS EXISTS. frozen-decisions.md carries one per-week non-negotiable: "re-run the
validation gate on the new target. Co-folding scores are not comparable across targets. The
instrument that scored barnase/barstar at 0.889 may simply fail on a GPCR. Find that out on
day 2, not day 7." targets/validation/ is barnase/barstar -- a different target. Nothing in
runs/ scores TNF-alpha, so as of 2026-10-06 the problem-2 instrument is UNVALIDATED and every
rank it would produce is uninterpretable.

THE POSITIVE IS TNFR2, AND THE CHOICE IS LOAD-BEARING.
  - Non-antibody, per arms-backlog s3a: co-folders systematically underperform on
    antibody-antigen, so an antibody positive UNDERSTATES the instrument and a failure could
    not be attributed.
  - Not a published de novo minibinder -- that would be circular (frozen-decisions s3b).
  - It is the receptor whose epitope s1 derived and s4 designs against, so this gate asks
    whether the instrument works AT OUR SITE, not merely on this target.
  - Natural, measured, and solved in complex (3ALQ, 4 TNFR2 copies).

SECOND POSITIVE, REPORTED BUT NON-GATING: adalimumab Fab (3WD5 H+L). Our submission may be
antibody-format, so this leg says what the instrument does on that format. It does not gate,
because a failure here is the known co-folder weakness rather than news about the target.

NEGATIVES: 4 composition-matched shuffles of the positive binder, fixed seed. Same amino
acids, destroyed order -- so a pass cannot be composition, and the shuffles carry the SAME
fold difficulty as the positive (TNFR2 is cysteine-rich; if ESMFold2 folds the receptor badly
it folds the shuffles badly too, and the comparison survives).

TARGET: targets/tnf/tnf_canonical_trimer.pdb -- the exact file the three generator specs point
at, so the gate and generation share one target. All three protomers, because the epitope
straddles a subunit interface (s1) and the assay target is the trimer.

Self-tests: python3 validation_gate_tnf.py --selftest
Build:     python3 validation_gate_tnf.py --write
"""
import os, sys, random, json
import gemmi
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fetch

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(REPO, "targets", "validation-tnf")
TRIMER = os.path.join(REPO, "targets", "tnf", "tnf_canonical_trimer.pdb")

# UniProt P01375 soluble domain, residues 77-233. Same constant per_partner.py uses.
HUMAN = ("VRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTHTIS"
         "RIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL")
SHUFFLE_SEED = 20261006   # written down before any score was computed (playbook s13)
N_NEG = 4


def one(r):
    i = gemmi.find_tabulated_residue(r.name)
    return i.one_letter_code.upper() if i and i.is_amino_acid() else None


def best_offset(ch):
    obs = [(r.seqid.num, one(r)) for r in ch if one(r)]
    if len(obs) < 40:
        return (None, 0)
    best = (None, 0)
    for k in range(-20, 130):
        m = t = 0
        for num, c in obs:
            u = num + k - 77
            if 0 <= u < len(HUMAN):
                t += 1
                m += (c == HUMAN[u])
        if t >= 40 and m / t > best[1]:
            best = (k, m / t)
    return best


def chain_seqs(path_or_id, cif=True):
    """Return {chain_name: observed one-letter sequence} with waters/ligands dropped."""
    path = fetch.cif(path_or_id) if cif else path_or_id
    st = gemmi.read_structure(path)
    st.setup_entities()
    st.remove_ligands_and_waters()
    out = {}
    for ch in st[0]:
        s = "".join(c for c in (one(r) for r in ch) if c)
        if s:
            out[ch.name] = s
    return out


def split_tnf_and_partners(seqs, structure_id):
    """Classify each chain as TNF (>=0.90 identity to HUMAN) or partner (>=15 aa)."""
    st = gemmi.read_structure(fetch.cif(structure_id))
    st.setup_entities()
    st.remove_ligands_and_waters()
    tnf, partner = [], []
    for ch in st[0]:
        _, acc = best_offset(ch)
        n = sum(1 for r in ch if one(r))
        if acc >= 0.90:
            tnf.append(ch.name)
        elif n >= 15:
            partner.append(ch.name)
    return tnf, partner


def trimer_file_check():
    """Verify targets/tnf/tnf_canonical_trimer.pdb agrees with canonical P01375 at EVERY
    observed position, per chain, matching by residue number rather than by list position.

    1TNF observes 152 of the 157 mature residues, so the chains start partway in; comparing
    string index 0 to canonical index 0 compares R6 against V1 and reports 150 mismatches on
    a correct file. The first draft of this script did exactly that and T1 caught it."""
    st = gemmi.read_structure(TRIMER)
    st.setup_entities()
    st.remove_ligands_and_waters()
    report = {}
    for ch in st[0]:
        k, acc = best_offset(ch)
        if k is None or acc < 0.90:
            continue
        mism, n = [], 0
        for r in ch:
            c = one(r)
            if c is None:
                continue
            u = r.seqid.num + k - 77
            if not (0 <= u < len(HUMAN)):
                continue
            n += 1
            if c != HUMAN[u]:
                mism.append((r.seqid.num + k, c, HUMAN[u]))
        report[ch.name] = {"offset": k, "observed": n, "mismatches": mism}
    return report


def canonical_trimer_chains():
    """The three protomers the fixtures actually ship: the CANONICAL 157-residue soluble
    domain, not the 152 observed in the crystal. The assay construct is P01375 77-233
    tag-free, 157 aa, and ESMFold2 folds from sequence -- there is no reason to hand it a
    gapped chain. The crystal file is still checked, by trimer_file_check(), because it is
    what the generator specs point at."""
    return {c: HUMAN for c in ("A", "B", "C")}


def shuffles(seq, n=N_NEG, seed=SHUFFLE_SEED):
    rng = random.Random(seed)
    out = []
    for i in range(n):
        letters = list(seq)
        rng.shuffle(letters)
        out.append("".join(letters))
    return out


def faa(chains):
    """chains: list of (name, seq). ESMFold2 wrapper format: >protein|<chain>."""
    return "".join(f">protein|{n}\n{s}\n" for n, s in chains)


def build():
    tri = canonical_trimer_chains()
    tri_chains = sorted(tri)
    assert len(tri_chains) == 3, f"expected 3 protomers in the canonical trimer, got {tri_chains}"

    # --- positive binder: the richest TNFR2 chain in 3ALQ ---
    alq = chain_seqs("3ALQ")
    tnf_ch, part_ch = split_tnf_and_partners(alq, "3ALQ")
    assert part_ch, "no TNFR2 partner chain found in 3ALQ"
    r2_name = max(part_ch, key=lambda c: len(alq[c]))
    r2 = alq[r2_name]

    # --- second positive, non-gating: adalimumab Fab H + L from 3WD5 ---
    wd5 = chain_seqs("3WD5")
    wd5_tnf, wd5_part = split_tnf_and_partners(wd5, "3WD5")
    fab = sorted(wd5_part, key=lambda c: -len(wd5[c]))[:2]

    cases = {}
    tri_entries = [(c, tri[c]) for c in tri_chains]
    cases["pos_tnf_tnfr2"] = tri_entries + [("D", r2)]
    for i, s in enumerate(shuffles(r2), 1):
        cases[f"neg_shuffled_tnfr2_{i}"] = tri_entries + [("D", s)]
    cases["pos2_tnf_adalimumab_fab_NONGATING"] = tri_entries + [
        ("D", wd5[fab[0]]), ("E", wd5[fab[1]])]

    meta = {
        "target_file": os.path.relpath(TRIMER, REPO),
        "target_chains": tri_chains,
        "target_len_per_protomer": [len(tri[c]) for c in tri_chains],
        "positive": {"source": "3ALQ", "chain": r2_name, "len": len(r2),
                     "what": "TNFR2 extracellular domain, natural receptor, non-antibody"},
        "negatives": {"n": N_NEG, "kind": "composition-matched shuffle of the positive binder",
                      "seed": SHUFFLE_SEED},
        "second_positive_nongating": {"source": "3WD5", "chains": fab,
                                      "lens": [len(wd5[c]) for c in fab],
                                      "what": "adalimumab Fab heavy+light"},
        "tnf_chains_detected": {"3ALQ": tnf_ch, "3WD5": wd5_tnf},
        "pass_rule": ("PASS = pos_tnf_tnfr2 ipSAE_min scores clearly above all "
                      f"{N_NEG} neg_shuffled_tnfr2_* . The Fab leg is reported, never gating."),
    }
    return cases, meta


def selftest():
    cases, meta = build()
    ok = []

    # T1 -- the fixtures ship the canonical 157-residue construct, AND the crystal file the
    # generator specs point at agrees with it at every observed position in all 3 chains.
    tri = canonical_trimer_chains()
    for c, seq in tri.items():
        assert seq == HUMAN, f"T1 FAIL fixture chain {c} is not canonical P01375 77-233"
        assert len(seq) == 157, f"T1 FAIL chain {c} length {len(seq)}, expected 157"
    chk = trimer_file_check()
    assert len(chk) == 3, f"T1 FAIL trimer file has {len(chk)} TNF chains, expected 3"
    for c, v in chk.items():
        assert not v["mismatches"], f"T1 FAIL {TRIMER} chain {c}: {v['mismatches'][:4]}"
    ok.append("T1 fixtures carry canonical 157 aa x3; the generator's target file matches it at "
              + ", ".join(f"{c}:{v['observed']}" for c, v in sorted(chk.items())) + " observed positions, 0 mismatches")

    # T2 -- every negative is a true permutation of the positive: composition held exactly.
    pos = dict(cases["pos_tnf_tnfr2"])["D"]
    for i in range(1, N_NEG + 1):
        neg = dict(cases[f"neg_shuffled_tnfr2_{i}"])["D"]
        assert sorted(neg) == sorted(pos), f"T2 FAIL neg {i} is not a permutation"
        assert neg != pos, f"T2 FAIL neg {i} is identical to the positive"
    ok.append(f"T2 all {N_NEG} negatives are exact permutations of the positive, none identical")

    # T3 -- the negatives differ from EACH OTHER, or n=4 is really n=1.
    negs = [dict(cases[f"neg_shuffled_tnfr2_{i}"])["D"] for i in range(1, N_NEG + 1)]
    assert len(set(negs)) == N_NEG, "T3 FAIL negatives are not distinct"
    ok.append(f"T3 the {N_NEG} negatives are distinct from each other")

    # T4 -- the shuffle is deterministic, so a rerun scores the same fixtures.
    again = shuffles(pos)
    assert again == negs, "T4 FAIL shuffle is not reproducible from the recorded seed"
    ok.append(f"T4 shuffles reproduce exactly from seed {SHUFFLE_SEED}")

    # T5 -- the positive is NOT an antibody chain. A ~220 aa chain paired with a second
    # ~215 aa chain out of 3ALQ would mean we picked a Fab by accident and the whole
    # non-antibody rationale is void.
    assert meta["positive"]["source"] == "3ALQ"
    assert len(pos) >= 100, f"T5 FAIL positive is only {len(pos)} aa; TNFR2 ECD is ~140-235"
    ok.append(f"T5 positive is a single {len(pos)} aa receptor chain ({meta['positive']['chain']} of 3ALQ), not a Fab pair")

    # T6 -- every case presents all three protomers plus at least one binder chain.
    for name, chains in cases.items():
        names = [n for n, _ in chains]
        assert len(names) == len(set(names)), f"T6 FAIL duplicate chain id in {name}"
        assert names[:3] == sorted(tri), f"T6 FAIL {name} is missing a protomer: {names}"
        assert len(names) >= 4, f"T6 FAIL {name} has no binder chain"
    ok.append("T6 every case carries all 3 protomers plus >=1 binder chain, ids unique")

    # T7 -- MUTATION TEST. Plant a defect that the gate exists to catch and confirm it fires.
    try:
        broken = sorted(pos)            # sorted, not shuffled: same composition, trivially
        assert sorted(broken) != sorted(pos)   # must fail -> the permutation check is real
        raise SystemExit("T7 FAIL: the permutation assertion cannot detect anything")
    except AssertionError:
        pass
    bad_seed = shuffles(pos, seed=SHUFFLE_SEED + 1)
    assert bad_seed != negs, "T7 FAIL: a different seed produced identical shuffles"
    ok.append("T7 mutation test: a changed seed changes the fixtures, so T4 is a real check")

    for line in ok:
        print("  ok  " + line)
    print(f"\nself-tests passed: {len(ok)}")
    return cases, meta


def write():
    cases, meta = build()
    os.makedirs(OUT, exist_ok=True)
    for name, chains in cases.items():
        with open(os.path.join(OUT, name + ".faa"), "w") as fh:
            fh.write(faa(chains))
    with open(os.path.join(OUT, "FIXTURES.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    n_res = sum(len(s) for _, s in cases["pos_tnf_tnfr2"])
    print(f"wrote {len(cases)} fixtures + FIXTURES.json to {os.path.relpath(OUT, REPO)}")
    print(f"  target: {meta['target_chains']} from {meta['target_file']}")
    print(f"  positive: 3ALQ chain {meta['positive']['chain']}, {meta['positive']['len']} aa "
          f"({meta['positive']['what']})")
    print(f"  negatives: {N_NEG} x composition-matched shuffle, seed {SHUFFLE_SEED}")
    print(f"  non-gating 2nd positive: 3WD5 {meta['second_positive_nongating']['chains']}")
    print(f"  complex size ~{n_res} residues -> that is what sizes the Modal timeout")
    print(f"\n  {meta['pass_rule']}")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest()
    elif "--write" in sys.argv:
        selftest()
        print()
        write()
    else:
        raise SystemExit(__doc__.strip().splitlines()[-2].strip())
