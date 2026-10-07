#!/usr/bin/env python3
"""Recompute the human->mouse epitope mapping from sequence and assert the JSON matches.

WHY THIS EXISTS. METHODS s2 claimed for days that human positional 68 is DELETED in mouse, and
built an argument on it: that a backbone difference made sidechain selection insufficient. It is
not deleted -- it is H->Y. The alignment has exactly ONE gap and it lands on human mature 71
(Ser), two residues upstream and not an epitope position. The nine mouse indices were right; the
ATTRIBUTION was wrong, and prose is where it hid, because no check recomputed it.

The lesson from s11 applies exactly: a mechanical check does not get tired or confident. This one
recomputes every mouse index, the gap location and every divergent-position label from the
sequences, two independent ways, and fails if the JSON or the prose disagrees.

Run: python3 bin/check_species_map.py
"""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "analysis", "02-tnf"))
from species_and_histidines import HUMAN, MOUSE, align

AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
       "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
       "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
JSON = os.path.join(ROOT, "analysis", "02-tnf", "epitope_conserved.json")
METHODS = os.path.join(ROOT, "submissions", "02-tnf-METHODS.md")


def observed(path, chain):
    d = {}
    with open(path) as f:
        for l in f:
            if l.startswith("ATOM") and l[21] == chain and l[12:16].strip() == "CA":
                d[int(l[22:26])] = AA3.get(l[17:20].strip(), "X")
    return "".join(d[k] for k in sorted(d))


def map_observed():
    """human positional -> (human res, mouse res or None, mouse positional or None)."""
    H = observed(os.path.join(ROOT, "targets", "tnf", "tnf_trimer_renum.pdb"), "A")
    M = observed(os.path.join(ROOT, "targets", "tnf", "tnf_mouse_trimer_renum.pdb"), "A")
    A, B = align(H, M)
    hp = mp = 0
    out = {}
    for x, y in zip(A, B):
        if x != "-": hp += 1
        if y != "-": mp += 1
        if x != "-":
            out[hp] = (x, (y if y != "-" else None), (mp if y != "-" else None))
    return out


def gaps_full():
    """Mature human positions with no mouse partner, from the FULL mature sequences."""
    A, B = align(HUMAN, MOUSE)
    m, out = 0, []
    for x, y in zip(A, B):
        if x != "-":
            m += 1
            if y == "-":
                out.append((m, x))
    return out


def main():
    d = json.load(open(JSON))
    obs = map_observed()
    fails = []

    # 1. every mouse epitope index must recompute
    for h, want in zip(d["human_positional"], d["mouse_positional"]):
        got = obs[h][2]
        if got != want:
            fails.append(f"epitope: human positional {h} maps to mouse {got}, JSON says {want}")

    # 2. every divergent label must match the aligned residues, and "deleted" must be real
    for pos, label in d["dropped_divergent"].items():
        hres, mres, _ = obs[int(pos)]
        if label.endswith("deleted"):
            if mres is not None:
                fails.append(f"dropped_divergent[{pos}] says deleted, but mouse has {mres}")
            continue
        m = re.fullmatch(r"([A-Z])->([A-Z])", label)
        if not m:
            fails.append(f"dropped_divergent[{pos}] label {label!r} is not parseable")
        elif (m.group(1), m.group(2)) != (hres, mres):
            fails.append(f"dropped_divergent[{pos}] says {label}, sequences say "
                         f"{hres}->{mres or 'DELETED'}")

    # 3. the gap must be where the JSON says, and there must be exactly as many as claimed
    g = gaps_full()
    md = d.get("mouse_deletion")
    if md is None:
        fails.append("JSON carries no mouse_deletion block; the gap location is unrecorded")
    else:
        if len(g) != md["gaps_in_alignment"]:
            fails.append(f"alignment has {len(g)} gap(s), JSON claims {md['gaps_in_alignment']}")
        if g and (g[0][0] != md["human_mature"] or g[0][1] != md["human_residue"]):
            fails.append(f"gap is at mature {g[0][0]} ({g[0][1]}), JSON says "
                         f"{md['human_mature']} ({md['human_residue']})")
        if md["human_positional"] != md["human_mature"] - 5:
            fails.append("mouse_deletion positional/mature are inconsistent (positional = mature-5)")

    # 4. the gap must not be an epitope position -- that was the substance of the old error
    if g and (g[0][0] - 5) in d["human_positional"]:
        fails.append(f"the gap at positional {g[0][0]-5} IS an epitope position; s2 must say so")

    # 5. Region I, if recorded, must recompute too
    r1 = d.get("region1")
    if r1:
        for h, want in zip(r1["human_positional"], r1["mouse_positional"]):
            got = obs[h][2]
            if got != want:
                fails.append(f"region1: human positional {h} maps to mouse {got}, "
                             f"JSON says {want}")
        a = r1.get("anchor_candidate", {})
        if a and obs[a["human_positional"]][2] != a["mouse_positional"]:
            fails.append("region1 anchor_candidate mouse_positional does not recompute")
        if a and obs[a["human_positional"]][1] != a["residue_both"]:
            fails.append("region1 anchor is not the same residue in mouse")

    # 6. METHODS must not reassert the withdrawn claim
    t = open(METHODS).read()
    for banned in ("mouse has no corresponding residue", "deleted in mouse**"):
        if banned in t:
            fails.append(f"METHODS still asserts {banned!r}; that claim was withdrawn")

    if fails:
        print("FAIL: species map does not reconcile with the sequences")
        for f in fails:
            print("  -", f)
        raise SystemExit(1)

    print(f"  {len(d['human_positional'])} epitope indices recomputed from the observed chains")
    print(f"  {len(d['dropped_divergent'])} divergent labels match the aligned residues")
    print(f"  exactly {len(g)} gap, at human mature {g[0][0]} ({g[0][1]}) = positional "
          f"{g[0][0]-5}, not an epitope position")
    if r1:
        print(f"  region1 {len(r1['human_positional'])} positions + anchor recomputed")
    print("\nPASS: every mouse index, the gap location and every divergent label "
          "recompute from sequence.")


if __name__ == "__main__":
    main()
