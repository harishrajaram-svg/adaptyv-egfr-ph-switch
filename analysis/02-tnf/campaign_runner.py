#!/usr/bin/env python3
"""Run any SIpHAB benchmark campaign through the repacked filter. 81 more single-point variants.

s18 ran Pertuzumab: 5 of 9, against a trivial null of 6 of 9. n=9 with 3 hits cannot carry that
conclusion. This generalises the run to every campaign with a deposited complex.

THE MAPPING PROBLEM, AND THE RULE THAT KEEPS IT HONEST
    The benchmark gives Kabat CDR positions (H31, L55, H100A). Structures are deposited in
    whatever the depositor used. 3BE1 matches Kabat exactly, 8 of 8. 5TRU and 8J6F are
    SEQUENTIAL, and a single global offset recovers only 11/21 and 33/52.

    Sequential and Kabat differ by a CONSTANT WITHIN A CDR and by a different constant in the
    next one, because the insertions sit between them. So the offset is fitted PER CHAIN PER
    CDR, and -- this is the part that matters -- a position is accepted ONLY if the residue
    identity at the mapped location equals the wild-type residue the benchmark names.

    Positions that do not verify are DROPPED, never guessed. s18's H100A was the one inferred
    mapping in that run and it produced a false positive; it is not a coincidence worth
    repeating. The count of dropped positions is reported with every campaign, because a
    campaign scored on 6 of 10 variants is not the campaign.

    \U0001F534 PDBFixer CANNOT ADDRESS INSERTION CODES. Chains are renumbered sequentially before
    any mutation (s18: asking for residue 99 found the residue at 99B instead, and 99 was a
    labelled hit). The Kabat map is kept through the renumbering.

Usage: campaign_runner.py <campaign key> [...]
"""
import collections, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dddg_elec_repacked as dd
import fetch

AA = {'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY','H':'HIS',
      'I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP',
      'Y':'TYR','V':'VAL'}
BENCH = os.path.join(HERE, "external", "sipHAB", "benchmark_parsed.json")
CDR_BANDS = {"cdr1": (24, 40), "cdr2": (45, 70), "cdr3": (89, 115)}

# campaign -> (pdb, {benchmark side: structure chain}, antigen chains)
CAMPAIGNS = {
    "bH1 (Her2)":          ("3BE1", {"H": "H", "L": "L"}, "A"),
    "Ipilimumab (CTLA-4)": ("5TRU", {"H": "H", "L": "L"}, "C"),
    "Tocilizumab (IL-6R)": ("8J6F", {"H": "H", "L": "L"}, "I"),
    "Pertuzumab (Her2)":   ("1S78", {"H": "D", "L": "C"}, "A"),
}


def parse_pos(pos):
    side, rest = pos[0], pos[1:]
    return side, int("".join(c for c in rest if c.isdigit())), "".join(c for c in rest if c.isalpha())


def band_of(num):
    for b, (lo, hi) in CDR_BANDS.items():
        if lo <= num <= hi:
            return b
    return "fw"


def prepare(pid, chains, antigen):
    """Protein-only, wanted chains, renumbered 1..N. Returns (path, renum, idx_by_orig)."""
    import gemmi
    st = gemmi.read_structure(fetch.cif(pid))
    st.setup_entities(); st.remove_ligands_and_waters(); st.remove_alternative_conformations()
    while len(st) > 1:
        del st[1]
    keep = set(chains.values()) | set(antigen)
    for cn in [c.name for c in st[0]]:
        if cn not in keep:
            st[0].remove_chain(cn)
    orig = {(c.name, r.seqid.num, r.seqid.icode.strip()): r.name for c in st[0] for r in c}
    renum = {}
    for ch in st[0]:
        for i, r in enumerate(ch, 1):
            renum[(ch.name, r.seqid.num, r.seqid.icode.strip())] = i
    for ch in st[0]:
        for i, r in enumerate(ch, 1):
            r.seqid.num = i; r.seqid.icode = " "
    st.setup_entities()
    out = os.path.join(HERE, f"{pid}_clean.pdb")
    st.write_pdb(out)
    return out, renum, orig


def map_positions(rows, chains, orig):
    """Per chain per CDR, fit one offset; accept a position only if the residue identity matches."""
    by = collections.defaultdict(list)
    for pos, mut, hit in rows:
        side, num, ic = parse_pos(pos)
        by[(side, band_of(num))].append((pos, side, num, ic, mut.split("->")[0], hit))
    mapped, dropped, offsets = [], [], {}
    for key, items in sorted(by.items()):
        side = key[0]
        best = (-1, 0)
        for off in range(-15, 16):
            ok = sum(1 for _, s, n, ic, wt, _h in items
                     if orig.get((chains[s], n + off, ic)) == AA.get(wt))
            if ok > best[0]:
                best = (ok, off)
        ok, off = best
        offsets[key] = (off, ok, len(items))
        for pos, s, n, ic, wt, hit in items:
            if orig.get((chains[s], n + off, ic)) == AA.get(wt):
                mapped.append((pos, chains[s], n + off, ic, wt, hit))
            else:
                dropped.append(pos)
    return mapped, dropped, offsets


def build(path, ch, num, wt, label):
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    out = os.path.join(HERE, f"_var_{label}.pdb")
    fx = PDBFixer(filename=path); fx.missingResidues = {}
    fx.applyMutations([f"{AA[wt]}-{num}-HIS"], ch)
    fx.findMissingAtoms(); fx.addMissingAtoms()
    with open(out, "w") as fh:
        PDBFile.writeFile(fx.topology, fx.positions, fh, keepIds=True)
    return out


def run(name):
    bench = json.load(open(BENCH))[name]
    pid, chains, antigen = CAMPAIGNS[name]
    binder = "".join(sorted(set(chains.values())))
    print(f"\n{'='*94}\n### {name}   {pid}   binder {binder} / antigen {antigen}")
    print(f"    benchmark: {bench['n']} single-point variants, {bench['hits']} labelled hits")
    path, renum, orig = prepare(pid, chains, antigen)
    mapped, dropped, offsets = map_positions(bench["rows"], chains, orig)
    print("    per-CDR offsets fitted (verified by residue identity):")
    for (side, band), (off, ok, n) in sorted(offsets.items()):
        print(f"      {side} {band:<5} offset {off:+3d}   {ok}/{n} verified")
    if dropped:
        print(f"    \U0001F534 DROPPED {len(dropped)} position(s) that did not verify: {dropped}")
    print(f"    scoring {len(mapped)} of {bench['n']} variants, "
          f"{sum(1 for m in mapped if m[5])} of {bench['hits']} hits retained")

    res = {}
    r = dd.score(path, binder, antigen, legs_wanted=("rigid", "repacked_local"))
    print(f"\n    {'variant':<10}{'labelled':<10}{'rigid':>9}{'repacked med [min..max]':>30}  pred")
    print("    " + "-" * 76)
    print(f"    {'WT':<10}{'-':<10}{r['rigid']['dddG']:>+9.3f}"
          f"{r['repacked_local']['dddG']:>+11.3f}")
    for pos, ch, num, ic, wt, hit in mapped:
        vp = build(path, ch, renum[(ch, num, ic)], wt, f"{pid}_{pos}")
        try:
            rv = dd.score(vp, binder, antigen, legs_wanted=("rigid", "repacked_local"))
        except SystemExit as e:
            print(f"    {pos:<10}REFUSED: {str(e)[:50]}"); continue
        rl = rv["repacked_local"]
        pred = rl["dddG"] >= 0
        res[pos] = (hit, pred, rl["dddG"])
        print(f"    {pos:<10}{('HIT' if hit else '-'):<10}{rv['rigid']['dddG']:>+9.3f}"
              f"{rl['dddG']:>+11.3f} [{rl['dddG_min']:+.2f}..{rl['dddG_max']:+.2f}]"
              f"{'!' if rl['sign_unstable'] else ' '}  {'switch' if pred else 'non':<7}"
              f"{'ok' if pred == hit else 'WRONG'}")
        os.remove(vp)
    if res:
        n = len(res); hits = sum(1 for h, _, _ in res.values() if h)
        ok = sum(1 for h, p, _ in res.values() if h == p)
        tp = sum(1 for h, p, _ in res.values() if h and p)
        fn = sum(1 for h, p, _ in res.values() if h and not p)
        fp = sum(1 for h, p, _ in res.values() if not h and p)
        null = max(hits, n - hits)
        print("    " + "-" * 76)
        print(f"    {ok}/{n} correct;  trivial null {null}/{n};  TP {tp}  FN {fn}  FP {fp}")
    return res


if __name__ == "__main__":
    keys = sys.argv[1:] or ["bH1 (Her2)"]
    allres = {}
    for k in keys:
        allres[k] = run(k)
    if len(allres) > 1:
        flat = [(h, p) for r in allres.values() for h, p, _ in r.values()]
        n = len(flat); hits = sum(1 for h, _ in flat if h)
        ok = sum(1 for h, p in flat if h == p)
        print(f"\n{'='*94}\nPOOLED over {len(allres)} campaigns: {ok}/{n} correct; "
              f"trivial null {max(hits, n-hits)}/{n};  "
              f"TP {sum(1 for h,p in flat if h and p)}  "
              f"FN {sum(1 for h,p in flat if h and not p)}  "
              f"FP {sum(1 for h,p in flat if not h and p)}")
