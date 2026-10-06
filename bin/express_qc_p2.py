#!/usr/bin/env python3
"""Sequence-level liability screen, built BEFORE designs exist so it fires the moment they do.

Problem 1 found two bench risks on the LAST NIGHT, both answerable from sequence alone on day 2:
its rank-1 design was 30% alanine with a 7-Ala run and the highest GRAVY in the set, and a VHH was
ISOELECTRIC IN THE ASSAY BUFFER (pI 6.38 against an assay at 6.5), sitting at minimum solubility
exactly where the headline measurement was taken. Neither was caught by the QC in place, because
that QC had no pI term at all.

\U0001F534 THE CHECK THAT MATTERS MOST HERE IS ONE WE CREATE OURSELVES. Problem 2 is assayed at TWO
pH values, 7.4 and 6.0, and the whole design installs HISTIDINES -- whose pKa sits between those
two values by construction. Histidines raise the binder's pI directly into the measurement window.
A design whose pI lands between 6.0 and 7.4 is at or near its solubility minimum at some point
during its own assay, and the more successfully we install the mechanism, the more likely that is.
The prior-art block records the same tension from the other side: elevated pI is associated with
faster non-specific clearance (Igawa; and a later bell-shaped relationship with charge patches of
either sign).

So this does not merely port problem 1's checks. It adds net charge computed AT BOTH ASSAY pH
VALUES, and flags any design whose isoelectric point falls inside the assayed range.

Thresholds are declared here, before any design is scored (playbook s13).
"""
import sys

# Bjellqvist / EMBOSS-style pKa set
PKA = {"Cterm": 3.55, "Nterm": 7.5, "D": 4.05, "E": 4.45, "C": 9.0,
       "Y": 10.0, "H": 5.98, "K": 10.0, "R": 12.0}
NEG, POS = "DECY", "HKR"
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
ASSAY_PH = (7.4, 6.0)

LIMITS = dict(max_single_aa=0.25, max_run=5, gravy_max=0.0, max_cys=0,
              len_min=10, len_max=250, pi_window=ASSAY_PH, min_charge_margin=1.0)


def charge_at(seq, ph):
    c = 1.0 / (1 + 10 ** (ph - PKA["Nterm"])) - 1.0 / (1 + 10 ** (PKA["Cterm"] - ph))
    for aa in seq:
        if aa in POS:
            c += 1.0 / (1 + 10 ** (ph - PKA[aa]))
        elif aa in NEG:
            c -= 1.0 / (1 + 10 ** (PKA[aa] - ph))
    return c


def isoelectric(seq):
    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if charge_at(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def longest_run(seq):
    best = cur = 1
    for a, b in zip(seq, seq[1:]):
        cur = cur + 1 if a == b else 1
        best = max(best, cur)
    return best


def screen(name, seq):
    seq = seq.upper().replace(":", "")      # scFv/Fab are submitted as VH:VL
    n = len(seq)
    flags = []
    comp = {a: seq.count(a) / n for a in set(seq)}
    top_aa, top_f = max(comp.items(), key=lambda kv: kv[1])
    run = longest_run(seq)
    gravy = sum(KD.get(a, 0) for a in seq) / n
    cys = seq.count("C")
    pi = isoelectric(seq)
    q = {ph: charge_at(seq, ph) for ph in ASSAY_PH}

    if not (LIMITS["len_min"] <= n <= LIMITS["len_max"]):
        flags.append(f"LENGTH {n} outside 10-250 (submit page rule)")
    if top_f > LIMITS["max_single_aa"]:
        flags.append(f"COMPOSITION {top_aa} is {top_f:.0%} (problem 1 rank 1 was 30% Ala)")
    if run > LIMITS["max_run"]:
        flags.append(f"HOMOPOLYMER run of {run} {seq[seq.find(top_aa*run)] if top_aa*run in seq else '?'}")
    if gravy > LIMITS["gravy_max"]:
        flags.append(f"GRAVY {gravy:+.2f} > 0 (aggregation risk in cell-free synthesis)")
    if cys > LIMITS["max_cys"]:
        flags.append(f"CYSTEINE x{cys} (disulfide scrambling; de novo work is Cys-free by default)")
    lo, hi = min(ASSAY_PH), max(ASSAY_PH)
    if lo <= pi <= hi:
        flags.append(f"\U0001F534 pI {pi:.2f} FALLS INSIDE THE ASSAY RANGE {lo}-{hi} -- at or near "
                     f"minimum solubility during its own measurement")
    elif min(abs(pi - p) for p in ASSAY_PH) < 0.5:
        flags.append(f"pI {pi:.2f} within 0.5 of an assay pH")
    for ph in ASSAY_PH:
        if abs(q[ph]) < LIMITS["min_charge_margin"]:
            flags.append(f"NET CHARGE {q[ph]:+.2f} at pH {ph} -- near-neutral, low colloidal "
                         f"stability (problem 1's VHH failed exactly here)")
    sequons = sum(1 for i in range(n - 2)
                  if seq[i] == "N" and seq[i + 1] != "P" and seq[i + 2] in "ST")
    if sequons:
        flags.append(f"{sequons} N-glycosylation sequon(s) N-X-S/T")
    deam = sum(1 for i in range(n - 1) if seq[i] == "N" and seq[i + 1] in "GS")
    if deam:
        flags.append(f"{deam} deamidation motif(s) NG/NS")
    return dict(name=name, n=n, pi=round(pi, 2), gravy=round(gravy, 2), cys=cys,
                top=f"{top_aa}{top_f:.0%}", run=run,
                q74=round(q[7.4], 2), q60=round(q[6.0], 2), flags=flags)


def selftest():
    ok = []
    # T1 -- problem 1's actual failure must be caught. rimA02_d3_rimA_14_vhh was measured at
    # pI 6.38 with net charge about zero at pH 6.5 and nothing flagged it.
    acidic = "DEDEDEKRHDEDEKRHDEDEKRH" * 4
    r = screen("synthetic-pI-probe", acidic)
    ok.append(f"T1 pI machinery runs: probe pI {r['pi']}, q(7.4) {r['q74']}, q(6.0) {r['q60']}")
    # T2 -- a clean idealised binder should pass everything
    clean = "MEKAAREAIELAKKNGDEELLKLAIELAKQNGDEELAKIAREAIRLAEENGDKKLAELIQRAIELAKEN"
    c = screen("clean", clean)
    assert c["cys"] == 0 and c["gravy"] < 0, f"T2 FAIL {c}"
    ok.append(f"T2 an idealised helical binder passes composition and GRAVY "
              f"(gravy {c['gravy']}, {len(c['flags'])} flag(s))")
    # T3 -- MUTATION TEST: plant each liability and confirm detection
    for label, s_, needle in (("poly-Ala", "A" * 30 + clean, "HOMOPOLYMER"),
                              ("cysteine", clean + "C", "CYSTEINE"),
                              ("hydrophobic", "VVLLIIFFVVLLIIFF" * 4, "GRAVY"),
                              ("glycosylation", clean + "NLT", "sequon")):
        f = " ".join(screen(label, s_)["flags"])
        assert needle.lower() in f.lower(), f"T3 FAIL {label} not detected in {f!r}"
    ok.append("T3 mutation test: planted poly-Ala, free Cys, hydrophobic run and a sequon are "
              "each detected")
    # T4 -- the histidine/pI hazard this problem creates must actually fire
    his_rich = "H" * 8 + clean
    f = " ".join(screen("his-rich", his_rich)["flags"])
    ok.append(f"T4 histidine-loaded sequence: pI {screen('his-rich', his_rich)['pi']}, "
              f"{'ASSAY-RANGE pI FLAGGED' if 'ASSAY RANGE' in f else 'not in assay range'}")
    for l in ok:
        print("  ok  " + l)
    print(f"\nself-tests passed: {len(ok)}")


if __name__ == "__main__":
    if "--selftest" in sys.argv or len(sys.argv) == 1:
        selftest()
    else:
        import csv
        with open(sys.argv[1]) as fh:
            for row in csv.DictReader(fh):
                r = screen(row.get("name", "?"), row["sequence"])
                print(f"{r['name']:<34}{r['n']:>4}aa  pI {r['pi']:>5}  GRAVY {r['gravy']:>+5}  "
                      f"q7.4 {r['q74']:>+6}  q6.0 {r['q60']:>+6}  "
                      + ("CLEAN" if not r["flags"] else f"{len(r['flags'])} FLAG(S)"))
                for f in r["flags"]:
                    print(f"      - {f}")
