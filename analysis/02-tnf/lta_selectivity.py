#!/usr/bin/env python3
"""Would a binder at this epitope also hit LYMPHOTOXIN-ALPHA? The selectivity axis nobody checked.

WHY THIS MATTERS AND WHY IT IS NOT IN THE BRIEF
    Lymphotoxin-alpha (LT-alpha / TNF-beta, UniProt P01374) shares TNF-alpha's jellyroll fold,
    trimerises the same way, and binds the SAME TWO RECEPTORS -- TNFR1 and TNFR2 -- at the SAME
    groove between adjacent protomers. That is the groove s1 derived as our epitope and s4
    designs into.

    The competition scores binding to human TNF-alpha, mouse TNF-alpha, and a pH ratio. It does
    NOT score selectivity, so cross-reactivity with LT-alpha costs nothing on the scoreboard. It
    costs something in a methods document that a reviewer reads, and it is the project's declared
    thesis: selectivity as a DESIGN objective rather than a post-hoc observation.

    Nothing in this project has looked at it. Thirteen sections of analysis, and the one protein
    most likely to be confused with the target has never been named.

WHAT THIS COMPUTES
    LT-alpha's residue at each of the 21 consensus receptor-epitope positions, by alignment to
    TNF-alpha, with the two design anchors called out. Free: one UniProt fetch and an alignment
    that already exists in species_and_histidines.py.

    Conserved anchors  -> a binder built on them is PREDICTED to cross-react, and that is a
                          designable constraint, not a disclaimer.
    Divergent anchors  -> selectivity comes free from the mechanism, which is a claim worth
                          making in the write-up because it was not designed for.
"""
import os, sys, urllib.request
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from species_and_histidines import HUMAN, MOUSE, UNIPROT_START, align

LTA_ACC = "P01374"
CACHE = os.path.join(HERE, "structures", f"{LTA_ACC}.fasta")

# s1's consensus receptor epitope: present in all 4 TNFR2 copies AND all 6 TNFR1 copies.
CONSENSUS = [96, 97, 107, 108, 109, 149, 151, 153, 161, 162, 163, 166, 167, 173,
             189, 191, 219, 220, 221, 222, 225]
ANCHORS = {108: "R108 - secondary anchor (s6d: reachable essentially only via NH2)",
           166: "K166 - PRIMARY anchor (s6d 17.2% reach, s9 the cation that separates)"}
MOUSE_DIVERGENT = {96, 107, 149, 161, 173}      # s2's five core residues that differ in mouse


def lta_sequence():
    if not os.path.exists(CACHE):
        url = f"https://rest.uniprot.org/uniprotkb/{LTA_ACC}.fasta"
        sys.stderr.write(f"fetching {url}\n")
        urllib.request.urlretrieve(url, CACHE)
    if os.path.getsize(CACHE) == 0:
        os.remove(CACHE)
        raise SystemExit(f"REFUSING: {LTA_ACC} fetch returned an empty file; removed the stub "
                         f"rather than letting a later read see zero residues")
    lines = open(CACHE).read().splitlines()
    return "".join(l.strip() for l in lines if not l.startswith(">"))


def main():
    lta = lta_sequence()
    A, B = align(HUMAN, lta)
    ident = sum(1 for a, b in zip(A, B) if a == b and a != "-")
    alen = sum(1 for a, b in zip(A, B) if a != "-" and b != "-")
    print(f"TNF-alpha (P01375 {UNIPROT_START}-233, {len(HUMAN)} aa) vs LT-alpha "
          f"({LTA_ACC}, {len(lta)} aa full chain)")
    print(f"  global alignment identity: {ident}/{alen} = {ident / alen:.1%}\n")

    # map TNF uniprot position -> aligned LT-alpha residue
    pos, m = UNIPROT_START - 1, {}
    for a, b in zip(A, B):
        if a != "-":
            pos += 1
            m[pos] = b
    mouse_m = {}
    MA, MB = align(HUMAN, MOUSE)
    pos = UNIPROT_START - 1
    for a, b in zip(MA, MB):
        if a != "-":
            pos += 1
            mouse_m[pos] = b

    print(f"{'pos':>5}  {'human':^6}{'mouse':^7}{'LT-a':^6}  verdict")
    print("-" * 74)
    same = diff = 0
    for u in CONSENSUS:
        h = HUMAN[u - UNIPROT_START]
        mo = mouse_m.get(u, "?")
        l = m.get(u, "-")
        ok = (l == h)
        same += ok
        diff += not ok
        tag = "conserved" if ok else f"DIVERGENT ({h}->{l})"
        extra = ""
        if u in ANCHORS:
            extra = "   <<< " + ANCHORS[u]
        elif u in MOUSE_DIVERGENT:
            extra = "   (also differs in mouse)"
        print(f"{u:>5}  {h:^6}{mo:^7}{l:^6}  {tag}{extra}")

    print("-" * 74)
    print(f"{same} of {len(CONSENSUS)} consensus epitope residues are SHARED with LT-alpha "
          f"({same / len(CONSENSUS):.0%}); {diff} differ.\n")

    print("THE TWO ANCHORS:")
    for u, why in sorted(ANCHORS.items()):
        h, l = HUMAN[u - UNIPROT_START], m.get(u, "-")
        state = "CONSERVED -- cross-reaction predicted" if l == h else \
                f"DIVERGENT ({h}->{l}) -- selectivity comes free here"
        print(f"  {why}\n      LT-alpha has {l!r}: {state}")

    cons_anchor = sum(1 for u in ANCHORS if m.get(u) == HUMAN[u - UNIPROT_START])
    print()
    if cons_anchor == len(ANCHORS):
        print("\U0001F534 BOTH anchors are conserved in LT-alpha. A binder whose mechanism IS the")
        print("   contact to those two cations has no structural reason to prefer TNF-alpha, so")
        print("   cross-reactivity should be the DEFAULT expectation and stated as such. The")
        print("   selectivity, if any, would have to come from the rest of the footprint.")
    elif cons_anchor == 0:
        print("✅ NEITHER anchor is conserved. Selectivity falls out of the mechanism itself,")
        print("   which is a stronger claim than a designed-in filter and was not designed for.")
    else:
        print(f"\U0001F7E0 {cons_anchor} of {len(ANCHORS)} anchors conserved. Selectivity is partial "
              f"and rests on the divergent one.")
    print("\n⚠️  Sequence-level only. Shared residues at shared positions make cross-reaction")
    print("   PLAUSIBLE, not certain -- the backbone geometry of the groove differs between the")
    print("   two proteins and nothing here measures that. A co-fold against LT-alpha would.")


if __name__ == "__main__":
    main()
