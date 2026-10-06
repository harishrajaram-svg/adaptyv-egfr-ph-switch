#!/usr/bin/env python3
"""Has ANY published molecule ever met this specification?

s1 computed what ratio the spec demands at a given affinity. It never asked the other half:
does anything in the measured record actually sit inside that region?

THE SPEC, objective 1: "binds human TNF-alpha at pH 7.4 and shows NO DETECTABLE BINDING at
pH 6.0". With a quantification ceiling of 10 uM (analysis/02-tnf/linkage_ceiling_p2.py), "no
detectable binding" means KD(6.0) must exceed that ceiling. So

    required ratio  >  CEILING / KD(7.4)

A tighter binder at 7.4 needs a LARGER ratio to disappear at 6.0. That is the weak-binder paradox
(playbook s2) stated as a feasibility question rather than a design preference.

\U0001F534 WHAT IS NEW HERE is the second column: every pH-switch molecule in this project's
prior-art record, placed on the same axes. Not "what would we need" but "what has anyone done".
"""
CEILING_nM = 10_000.0          # 10 uM quantification ceiling

# (name, KD at the ON pH in nM, measured ratio, what KIND of ratio, pH pair, source)
RECORD = [
    ("Chugai anti-myostatin",   0.285, 119.0, "KD",       "7.4 -> 5.8",
     "10.1080/19420862.2022.2068213 -- best VERIFIED KD ratio in the record"),
    ("Satralizumab (marketed)", None,   10.0,  "KD",       "7.4 -> 6.0",
     "EMA EPAR, flagged UNVERIFIED in prior-art"),
    ("adalimumab WT",           0.0046, 9.0,   "off-rate", "7.4 -> 6.0",
     "Schroter 2015, 10.4161/19420862.2014.985993"),
    ("Schroter PSV#1",          0.0463, 231.0, "off-rate", "7.4 -> 6.0", "Schroter 2015"),
    ("Schroter PSV#2",          0.0773, 785.0, "off-rate", "7.4 -> 6.0", "Schroter 2015"),
    ("Schroter PSV#3",          0.112,  505.0, "off-rate", "7.4 -> 6.0", "Schroter 2015"),
    ("Adafre AF-M2637",         1.1,    15.0,  "EC50",     "7.4 -> 6.0",
     "Watkins 2022, PMID 35896334"),
    ("Adafre AF-M2631",         None,   30.0,  "EC50",     "7.4 -> 6.0", "Watkins 2022"),
    ("best computational, any", None,   2.0,   "KD",       "7.4 -> 6.0",
     "prior-art: best computational result anywhere in this window"),
]
OURS = [("s4 pessimistic, Broo 1 cation/site", 4.5), ("s4 optimistic, 2.0-unit drop", 51.6)]


def need(kd_nM):
    return CEILING_nM / kd_nM


def main():
    print(__doc__.strip().splitlines()[0] + "\n")
    print(f"Ceiling {CEILING_nM/1000:.0f} uM. Required ratio = ceiling / KD(7.4).\n")
    print(f"{'KD at 7.4':>12}{'ratio needed to vanish at 6.0':>32}")
    print("-" * 46)
    for kd in (0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 500.0, 1000.0, 5000.0):
        print(f"{kd:>10.3f} nM{need(kd):>30,.0f}x")

    print(f"\n{'molecule':<30}{'KD@7.4':>10}{'ratio':>9}{'kind':>10}{'needed':>12}  verdict")
    print("-" * 92)
    any_pass = False
    for name, kd, ratio, kind, pair, src in RECORD:
        if kd is None:
            print(f"{name:<30}{'n/a':>10}{ratio:>8.0f}x{kind:>10}{'--':>12}  "
                  f"KD not in hand; cannot place")
            continue
        req = need(kd)
        ok = ratio >= req
        any_pass |= ok
        print(f"{name:<30}{kd:>10.4f}{ratio:>8.0f}x{kind:>10}{req:>11,.0f}x  "
              f"{'MEETS SPEC' if ok else 'fails by ' + f'{req/ratio:,.0f}x'}")

    print("\n" + "=" * 92)
    if not any_pass:
        print("\U0001F534 NOTHING IN THE RECORD MEETS THE SPEC.")
        best = min((need(kd) / r, n, kd, r) for n, kd, r, k, p, s in RECORD if kd)
        print(f"   Closest is {best[1]}: it would need to be {best[0]:,.0f}x better.")
    print("""
WHY, in one line: the molecules people actually build are PICOMOLAR at the ON pH, and a picomolar
binder needs a four-to-seven-log ratio to cross a 10 uM ceiling. Nobody gets near that. The
record's best VERIFIED KD ratio is 119x; the spec wants 35,000x of a 285 pM binder.

THE ESCAPE, and it is the whole design consequence: bind WEAKLY at 7.4. At KD 500 nM the spec
needs 20x; at 1 uM it needs 10x. Both sit INSIDE the published 5-30x band. Our own s4 band of
4.5-51.6x clears 1 uM comfortably and 500 nM at its optimistic end.

So objective 1 is reachable -- but ONLY in a band where objective 3 (affinity) is deliberately
sacrificed. That is not a tension we introduced; it is in the arithmetic of the spec. It is also
exactly what s4 already chose (200 nM - 1 uM), which this now justifies from the opposite
direction: not "what can our mechanism deliver" but "what does the spec permit at all".""")
    print("""
⚠️  THREE QUANTITIES, NOT ONE. The ratios above are KD, off-rate and EC50 shifts and they are
   NOT interchangeable. Only the KD rows are strictly comparable to the requirement, which is a
   KD question. The Schroter and Adafre rows are placed for scale and are marked.

⚠️  The 10 uM ceiling is OUR assumption about the assay. Organiser question 2 asks what "no
   detectable binding" actually means. If it means "no response above buffer" at a lower
   concentration, every number above moves. That question is now the single highest-value
   unknown on this problem.""")


if __name__ == "__main__":
    main()
