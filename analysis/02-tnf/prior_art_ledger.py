"""Prior-art exclusion ledger for problem 2, keyed on BINDER SEQUENCE.

NOT the same thing as bin/exclusion_ledger.py. That one is problem 1's ledger of OUR OWN
designs, classified by why we cut them (playbook s38: "cut before measurement" vs "measured
and rejected"). This one is the opposite direction: every KNOWN TNF-alpha binder, so that
(a) no design we generate is an unwitting rediscovery, and (b) there is a local novelty
reference before the organisers' checker runs at upload. s5 of the challenge file says the
21-entry RCSB survey "is a starting point, not the novelty gate" -- this turns the starting
point into sequences, which is the only form a novelty check can actually use.

Problem 1 taught the specific lesson: our gate cleared 18 of 18 and Proteinbase rejected two
(playbook s34). A local ledger does not replace their checker. It catches the rediscoveries
their checker would catch AT UPLOAD, when there may be no retry left.

ENTRIES
  structure-derived  binder (non-TNF) chains from solved human TNF-alpha complexes, pulled
                     from RCSB. Certain, and reproducible in a clean clone via fetch.py.
  rebuilt            Schroter et al.'s PSV#1/2/3 adalimumab histidine variants, as rebuilt
                     for the s6b known-answer control. Prior art on our switch direction,
                     on our target, with measured ratios.
  publication-only   named, sequences NOT retrieved. Recorded as a GAP with the reason, never
                     as an empty row that reads like a clean check (s10: a zero is not an
                     absence).

Self-tests: validation via --selftest. Build: --write.
"""
import os, sys, json
import gemmi
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fetch

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(HERE, "prior_art_ledger.json")
TSV = os.path.join(HERE, "prior_art_ledger.tsv")

HUMAN = ("VRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTHTIS"
         "RIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL")

# The 12 with contacts computed in s5, plus the title-only group, plus one s5 missed.
STRUCTURES = {
    "3WD5": ("adalimumab Fab", "approved antibody"),
    "4G3Y": ("infliximab Fab", "approved antibody"),
    "5WUX": ("certolizumab Fab", "approved antibody"),
    "5YOY": ("golimumab Fab", "approved antibody"),
    "21TW": ("ozoralizumab / TNF30", "approved antibody, VHH-based"),
    "21TV": ("ozoralizumab-HSA", "approved antibody, VHH-based"),
    "8Z8M": ("TNF30 VHH", "nanobody"),
    "5M2I": ("llama VHH1", "nanobody"),
    "5M2J": ("llama VHH2", "nanobody"),
    "5M2M": ("llama VHH3", "nanobody"),
    "9BN7": ("VNAR C4", "shark VNAR"),
    "9DJW": ("VNAR D1", "shark VNAR -- 30-chain lattice, chain split unverified"),
    "3ALQ": ("TNFR2 ectodomain", "natural receptor"),
    "8ZUI": ("TNFR1 ectodomain", "natural receptor"),
    "7KPB": ("TNFR1 ectodomain", "natural receptor"),
    "3IT8": ("poxvirus TNF-binding protein", "viral decoy receptor"),
    "7TA3": ("synthetic alpha/beta-peptide", "peptidomimetic"),
    "7TA6": ("synthetic alpha/beta-peptide", "peptidomimetic"),
    "4TWT": ("bicyclic peptide M21", "parses as a LIGAND, not a polymer chain -- expect 0 chains"),
    # s5 does not list this one. Murtaugh et al., Protein Sci 20(9):1619-1631 (2011),
    # doi 10.1002/pro.696 -- combinatorial histidine scanning in a VHH, OUR format and OUR
    # mechanism, with a 1.75 A structure of the 5-histidine variant. The challenge file calls
    # it "highest value per unit effort" among leads not yet followed; it is also prior art.
    # \U0001F534 3QSK REMOVED 2026-10-06, same day it was added. It is NOT a TNF complex: its own
    # title is "5 Histidine Variant of the anti-RNase A VHH in Complex with RNAse A", and its
    # two chains are Ribonuclease pancreatic and an anti-RNase A camelid VHH. The challenge
    # file describes Murtaugh et al. as "combinatorial histidine scanning in a VHH, OUR FORMAT"
    # -- it never said our target, and this script read "our format" as "our target".
    #
    # The cost of leaving it in: both chains fail the >=0.90 TNF test, so BOTH were filed as
    # TNF-binder prior art. RNase A is pure noise, and the anti-RNase VHH would raise a false
    # rediscovery alarm against any nanobody sharing a common camelid framework -- which the
    # challenge FAQ explicitly permits. T7 below now makes this class of error impossible.
    #
    # Murtaugh remains valuable and is NOT lost: it is a real backbone carrying five installed
    # histidines with a 1.75 A structure, which is the best available fixture for asking whether
    # a BACKBONE holds a histidine in place where our free tripeptide probe does not (s10's
    # largest open confound). Different target, same question. Tracked in the challenge file,
    # not here.
}

REBUILT = {
    # analysis/02-tnf/structures/3WD5_PSV{1,2,3}.pdb, built for the s6b control.
    "3WD5_PSV1": ("Schroter PSV#1", "adalimumab + 5 His, measured kd ratio 231x"),
    "3WD5_PSV2": ("Schroter PSV#2", "adalimumab + 5 His, measured kd ratio 785x"),
    "3WD5_PSV3": ("Schroter PSV#3", "adalimumab + 3 His, measured kd ratio 505x"),
    # Added 2026-10-06 (s17). Adafre Biosciences, Watkins & Watkins, J Immunol 209(4):829-839
    # (2022), PMID 35896334, doi 10.4049/jimmunol.2101180. pH-sensitive monovalent adalimumab
    # variants -- OUR target, OUR direction, PATENT PENDING per the paper's conflict statement.
    # They post-date the 21-entry RCSB survey by seven years and nothing in this project had
    # named them until the SIpHAB reference list surfaced them.
    "3WD5_AFM2637": ("Adafre AF-M2637",
                     "adalimumab L chain Q89H/R90H/N92H, monovalent; EC50 15x higher after a "
                     "pH 6.0 wash. Patent pending"),
    "3WD5_AFM2631": ("Adafre AF-M2631",
                     "adalimumab H chain L102H -- ONE histidine, EC50 30x. The single-His "
                     "positive our filter called a non-switch (s17). Patent pending"),
}

# Named, sequences not in hand. Each carries WHY, so the gap is auditable.
PUBLICATION_ONLY = [
    {"name": "Ahn et al. TNF-alpha binder designs",
     "n_designs": 72, "n_ph_sensitive": 3,
     "source": "bioRxiv 10.1101/2025.09.29.678932 (2025-09-29), Baker lab",
     "why_a_gap": "no code released; sequences would be in a supplementary table that has not "
                  "been retrieved. This is prior art on our EXACT target AND our exact "
                  "mechanism (interface His against a cation), published after the 21-entry "
                  "RCSB survey, so it is the single most consequential gap in this ledger.",
     "blocks": "a design of ours could duplicate one of the 72 and we could not tell."},
    {"name": "Merck KGaA adalimumab pH-variant patent",
     "source": "US10183994B2 / WO2016000813A1",
     "why_a_gap": "claims adalimumab variants by kdis ratio >=5/10/15x over wild type. Claim "
                  "scope, not a sequence set; the PSV rows above are its published examples.",
     "blocks": "nothing for a de novo entry -- recorded as prior art to cite, not to avoid."},
]


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


def binder_chains(path):
    """Every chain in `path` that is NOT TNF-alpha and is >= 15 aa.

    The >=15 floor and the >=0.90 TNF-identity test are lifted from per_partner.py so the
    two agree on what counts as a partner chain; divergent definitions of 'binder' between
    two scripts is how a count stops meaning anything."""
    st = gemmi.read_structure(path)
    st.setup_entities()
    st.remove_ligands_and_waters()
    out, tnf = {}, []
    for ch in st[0]:
        _, acc = best_offset(ch)
        seq = "".join(c for c in (one(r) for r in ch) if c)
        if acc >= 0.90:
            tnf.append(ch.name)
        elif len(seq) >= 15:
            out[ch.name] = seq
    return out, tnf


def has_tnf(path):
    """Does this structure actually contain a TNF-alpha chain?

    Added after 3QSK -- an anti-RNase A VHH complex -- was filed as TNF prior art because the
    source note said "our format" and this script read it as "our target". Every chain failed
    the TNF test, so every chain was classified as a binder. A ledger that cannot tell whether
    the target is in the file will happily record binders to something else."""
    _, tnf = binder_chains(path)
    return bool(tnf)


def build():
    entries, problems = [], []

    for pid, (name, kind) in sorted(STRUCTURES.items()):
        try:
            path = fetch.cif(pid)
        except SystemExit as e:
            problems.append({"id": pid, "stage": "fetch", "detail": str(e)})
            continue
        try:
            chains, tnf = binder_chains(path)
        except Exception as e:
            problems.append({"id": pid, "stage": "parse", "detail": repr(e)})
            continue
        if not tnf:
            problems.append({"id": pid, "stage": "NOT-A-TNF-COMPLEX",
                             "detail": f"no chain reads as TNF-alpha; chains "
                                       f"{sorted(chains)} would have been filed as TNF "
                                       f"binders. REFUSED."})
            continue
        if not chains:
            problems.append({"id": pid, "stage": "no-binder-chain",
                             "detail": f"TNF chains {tnf}, 0 partner chains >=15 aa"})
        for cname, seq in sorted(chains.items()):
            entries.append({"key": seq, "pdb": pid, "chain": cname, "name": name,
                            "kind": kind, "len": len(seq), "basis": "structure",
                            "tnf_chains": tnf})

    for stem, (name, kind) in sorted(REBUILT.items()):
        path = os.path.join(HERE, "structures", stem + ".pdb")
        if not os.path.exists(path):
            problems.append({"id": stem, "stage": "missing-rebuild", "detail": path})
            continue
        chains, tnf = binder_chains(path)
        for cname, seq in sorted(chains.items()):
            entries.append({"key": seq, "pdb": stem, "chain": cname, "name": name,
                            "kind": kind, "len": len(seq), "basis": "rebuilt",
                            "tnf_chains": tnf})

    by_seq = {}
    for e in entries:
        by_seq.setdefault(e["key"], []).append(e)

    meta = {
        "generated": "2026-10-06",
        "n_entries": len(entries),
        "n_distinct_sequences": len(by_seq),
        "n_structures_attempted": len(STRUCTURES) + len(REBUILT),
        "n_structures_with_no_binder_chain": sum(1 for p in problems
                                                 if p["stage"] == "no-binder-chain"),
        "problems": problems,
        "publication_only_GAPS": PUBLICATION_ONLY,
        "how_to_use": ("Keyed on binder sequence. Before the novelty gate runs, every generated "
                       "design is checked against these sequences. This is a LOCAL reference and "
                       "does NOT replace the organisers' checker at upload (playbook s34)."),
    }
    return entries, by_seq, meta


def identity(a, b):
    """Ungapped identity over the shorter sequence, at the best single offset.

    Deliberately crude and deliberately NOT the novelty gate. bin/novelty_gate.py runs
    foldseek with TMalign and a coverage clause; that is the instrument. This is a cheap
    rediscovery alarm for use during generation, when the question is only 'have we just
    regenerated adalimumab'. A real near-duplicate reads very high here; anything subtle is
    the gate's business, not this function's."""
    if not a or not b:
        return 0.0
    if len(a) > len(b):
        a, b = b, a
    best = 0.0
    for off in range(0, len(b) - len(a) + 1):
        m = sum(1 for i, c in enumerate(a) if b[off + i] == c)
        best = max(best, m / len(a))
    return best


def check(seq, by_seq=None, top=3):
    """Report the closest prior-art sequences to `seq`."""
    if by_seq is None:
        _, by_seq, _ = build()
    scored = sorted(((identity(seq, k), k) for k in by_seq), reverse=True)
    out = []
    for frac, k in scored[:top]:
        e = by_seq[k][0]
        out.append({"identity": round(frac, 4), "pdb": e["pdb"], "chain": e["chain"],
                    "name": e["name"], "len": e["len"]})
    return out


def selftest():
    entries, by_seq, meta = build()
    ok = []

    # T1 -- the ledger is not empty, and it is not trivially small. s10: a zero is not an
    # absence, and a ledger that silently fetched nothing would read as "no prior art".
    assert len(entries) >= 20, f"T1 FAIL only {len(entries)} entries"
    assert len(by_seq) >= 15, f"T1 FAIL only {len(by_seq)} distinct sequences"
    ok.append(f"T1 {len(entries)} binder chains, {len(by_seq)} distinct sequences, "
              f"from {meta['n_structures_attempted']} structures")

    # T2 -- no TNF chain leaked in as a 'binder'. If one did, every design would read ~100%
    # identical to prior art the moment it was checked against the target itself.
    for e in entries:
        assert identity(e["key"], HUMAN) < 0.90, \
            f"T2 FAIL {e['pdb']}:{e['chain']} is {identity(e['key'], HUMAN):.2f} identical to TNF"
    ok.append("T2 no entry is TNF-alpha itself (all < 0.90 identity to canonical)")

    # T3 -- KNOWN ANSWER. Adalimumab's own heavy chain must come back at identity 1.0 against
    # the ledger, and the PSV variants must come back HIGH BUT NOT 1.0 -- they are adalimumab
    # plus installed histidines. Both directions in one test.
    ada = [e for e in entries if e["pdb"] == "3WD5" and e["basis"] == "structure"]
    assert ada, "T3 FAIL no 3WD5 structure rows"
    hit = check(ada[0]["key"], by_seq, top=1)[0]
    assert hit["identity"] == 1.0, f"T3 FAIL adalimumab chain reads {hit['identity']}"
    psv = [e for e in entries if e["pdb"].startswith("3WD5_PSV")]
    if psv:
        longest = max(psv, key=lambda e: e["len"])
        best_struct = max(identity(longest["key"], e["key"]) for e in ada)
        assert 0.80 <= best_struct < 1.0, \
            f"T3 FAIL a PSV chain is {best_struct:.3f} against wild-type adalimumab"
        ok.append(f"T3 known answer: adalimumab reads 1.000 against itself; its PSV variant "
                  f"reads {best_struct:.3f} — high, and NOT identical, which is what "
                  f"installed histidines should look like")
    else:
        ok.append("T3 known answer: adalimumab reads 1.000 against itself (no PSV rows found)")

    # T4 -- a sequence that is NOT prior art must score low. Without this, T3 passing is
    # consistent with identity() returning 1.0 for everything.
    decoy = "M" + "EKAALRELAIKAVELAKQNGDEKLLELAIKLAEQNGDKELAEKIAREAIRLAEENGDKKLAELI"
    worst = check(decoy, by_seq, top=1)[0]["identity"]
    assert worst < 0.45, f"T4 FAIL a designed-looking decoy reads {worst:.3f} against prior art"
    ok.append(f"T4 negative control: an idealised helical decoy reads {worst:.3f} — "
              f"identity() is not saturated")

    # T5 -- MUTATION TEST. Plant a near-duplicate (one substitution into a real prior-art
    # sequence) and confirm the alarm fires at >=0.95. A ledger that cannot see a 1-residue
    # variant of adalimumab would not have caught problem 1's failure mode either.
    real = max(by_seq, key=len)
    mut = ("A" if real[0] != "A" else "G") + real[1:]
    fired = check(mut, by_seq, top=1)[0]["identity"]
    assert fired >= 0.95, f"T5 FAIL a 1-substitution near-duplicate reads only {fired:.3f}"
    ok.append(f"T5 mutation test: a 1-substitution near-duplicate of the longest entry reads "
              f"{fired:.3f}")

    # T7 -- EVERY entry comes from a structure that actually contains TNF-alpha. This is the
    # 3QSK guard: without it, a complex of something else entirely files both of its chains as
    # TNF-binder prior art, and the ledger's whole purpose -- "have we regenerated a known TNF
    # binder" -- is answered against the wrong molecule.
    for e in entries:
        assert e["tnf_chains"], f"T7 FAIL {e['pdb']}:{e['chain']} has no TNF chain in its source"
    ok.append(f"T7 all {len(entries)} entries come from structures containing a TNF-alpha chain; "
              f"{sum(1 for p in meta['problems'] if p['stage'] == 'NOT-A-TNF-COMPLEX')} "
              f"structure(s) refused for not being TNF complexes")

    # T6 -- the gaps are carried as gaps. If PUBLICATION_ONLY were ever emptied, the ledger
    # would read complete while missing 72 designs on our exact target.
    assert meta["publication_only_GAPS"], "T6 FAIL the publication-only gap list is empty"
    ahn = [g for g in meta["publication_only_GAPS"] if "Ahn" in g["name"]]
    assert ahn and ahn[0]["n_designs"] == 72, "T6 FAIL Ahn's 72 designs are not recorded as a gap"
    ok.append(f"T6 {len(meta['publication_only_GAPS'])} publication-only gaps carried explicitly, "
              f"Ahn's 72 designs among them")

    for line in ok:
        print("  ok  " + line)
    if meta["problems"]:
        print("\n  structures that yielded no usable binder chain (reported, not hidden):")
        for p in meta["problems"]:
            print(f"    {p['id']:<10} {p['stage']:<20} {p['detail'][:90]}")
    print(f"\nself-tests passed: {len(ok)}")
    return entries, by_seq, meta


def write():
    entries, by_seq, meta = build()
    json.dump({"meta": meta, "entries": entries}, open(OUT, "w"), indent=1)
    with open(TSV, "w") as fh:
        fh.write("pdb\tchain\tname\tkind\tlen\tbasis\tsequence\n")
        for e in sorted(entries, key=lambda x: (x["pdb"], x["chain"])):
            fh.write(f"{e['pdb']}\t{e['chain']}\t{e['name']}\t{e['kind']}\t{e['len']}\t"
                     f"{e['basis']}\t{e['key']}\n")
    print(f"wrote {os.path.relpath(OUT, REPO)} and {os.path.relpath(TSV, REPO)}")
    print(f"  {meta['n_entries']} binder chains, {meta['n_distinct_sequences']} distinct sequences")
    print(f"  {len(meta['publication_only_GAPS'])} publication-only GAPS carried, "
          f"Ahn's 72 designs the consequential one")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest()
    elif "--write" in sys.argv:
        selftest(); print(); write()
    elif "--check" in sys.argv:
        q = sys.argv[sys.argv.index("--check") + 1]
        for h in check(q):
            print(f"  {h['identity']:.4f}  {h['pdb']}:{h['chain']}  {h['name']} ({h['len']} aa)")
    else:
        raise SystemExit("usage: prior_art_ledger.py [--selftest|--write|--check SEQ]")
