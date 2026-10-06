#!/usr/bin/env python3
"""Invert the question: what do the 14 MEASURED switches have in common?

Every instrument in this project tests OUR hypothesis against the data. s19 ended that line --
across ~100 real molecules, a binder histidine accepting from a cross-interface cation has never
been observed in a confirmed switch, while the cations sit within 8 A in 66 of 82 cases.

So stop asking "does the data support our criterion" and ask the data what ITS criterion is.

We hold something almost nobody does: 82 single-point histidine variants on four antigens, each
with a measured hit / non-hit label (s16), and a deposited complex for every one. The design
question -- "which position should I mutate to histidine" -- is answerable from the PARENT
structure, before any mutation, which is exactly how a generator would have to use it.

FEATURES, all computed on the wild-type complex at the position about to be mutated:
    burial        heavy atoms within 10 A of the sidechain (the proxy third_site_census uses)
    d_antigen     min sidechain distance to any antigen heavy atom
    n_contact     antigen heavy atoms within 6 A
    d_anion       nearest antigen Asp/Glu carboxylate oxygen
    d_cation      nearest antigen Arg/Lys/His cationic nitrogen
    n_arom6       antigen aromatic ring atoms within 6 A
    wt_charge     is the residue being replaced charged
    chain, cdr

\U0001F534 THE HONEST FRAME. 14 positives against 68 negatives, four antigens, and one feature
will separate them by chance if enough are tried. Every split below is reported with the count
on BOTH sides and nothing is called a finding on a single campaign. s13: a threshold needs a bar
set before the numbers -- so the bar here is stated first: a feature is worth reporting only if
it separates on at least THREE of the four campaigns in the same direction.
"""
import collections, json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import campaign_runner as cr

ANION = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}
CATION = {("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"), ("LYS", "NZ"),
          ("HIS", "ND1"), ("HIS", "NE2")}
AROM = {"PHE", "TYR", "TRP", "HIS"}
BACKBONE = {"N", "CA", "C", "O"}
MIN_CAMPAIGNS = 3        # written down before any number was computed


def features(path, chains, antigen, ch, num):
    import gemmi
    st = gemmi.read_structure(path)
    st.setup_entities(); st.remove_ligands_and_waters()
    m = st[0]
    tgt, side, res = [], [], None
    for c in m:
        for r in c:
            for a in r:
                if a.element == gemmi.Element("H"):
                    continue
                p = np.array([a.pos.x, a.pos.y, a.pos.z])
                if c.name in antigen:
                    tgt.append((r.name, a.name, p))
                if c.name == ch and r.seqid.num == num:
                    res = r
                    if a.name not in BACKBONE:
                        side.append(p)
    if res is None or not side or not tgt:
        return None
    side = np.asarray(side)
    tp = np.asarray([p for _, _, p in tgt])
    d = np.linalg.norm(side[:, None, :] - tp[None, :, :], axis=2)
    dmin = float(d.min())

    allatoms = np.asarray([[a.pos.x, a.pos.y, a.pos.z] for c in m for r in c for a in r
                           if a.element != gemmi.Element("H")])
    burial = int((np.linalg.norm(allatoms - side.mean(axis=0), axis=1) <= 10.0).sum())

    def nearest(sel):
        pts = np.asarray([p for rn, an, p in tgt if (rn, an) in sel])
        if not len(pts):
            return 99.0
        return float(np.min(np.linalg.norm(side[:, None, :] - pts[None, :, :], axis=2)))

    arom = np.asarray([p for rn, an, p in tgt if rn in AROM and an not in BACKBONE])
    n_arom = int((np.linalg.norm(side[:, None, :] - arom[None, :, :], axis=2) <= 6.0).sum()) \
        if len(arom) else 0
    return dict(burial=burial, d_antigen=round(dmin, 2),
                n_contact=int((d <= 6.0).sum()),
                d_anion=round(nearest(ANION), 2), d_cation=round(nearest(CATION), 2),
                n_arom6=n_arom, wt=res.name)


def main():
    bench = json.load(open(cr.BENCH))
    rows = []
    for name in cr.CAMPAIGNS:
        pid, chains, antigen = cr.CAMPAIGNS[name]
        path, renum, orig = cr.prepare(pid, chains, antigen)
        mapped, dropped, _ = cr.map_positions(bench[name]["rows"], chains, orig)
        for pos, ch, num, ic, wt, hit in mapped:
            f = features(path, chains, antigen, ch, renum[(ch, num, ic)])
            if f is None:
                continue
            f.update(campaign=name, pos=pos, hit=bool(hit))
            rows.append(f)
        print(f"  {name:<24}{len([r for r in rows if r['campaign']==name]):>3} positions featurised")

    hits = [r for r in rows if r["hit"]]
    non = [r for r in rows if not r["hit"]]
    print(f"\n{len(rows)} positions: {len(hits)} measured switches, {len(non)} non-switches\n")

    print(f"{'feature':<12}{'hits median':>13}{'non median':>12}{'separates on':>14}  direction")
    print("-" * 74)
    keep = []
    for feat in ("burial", "d_antigen", "n_contact", "d_anion", "d_cation", "n_arom6"):
        mh = float(np.median([r[feat] for r in hits]))
        mn = float(np.median([r[feat] for r in non]))
        agree = 0
        for name in cr.CAMPAIGNS:
            h = [r[feat] for r in rows if r["campaign"] == name and r["hit"]]
            n = [r[feat] for r in rows if r["campaign"] == name and not r["hit"]]
            if not h or not n:
                continue
            if (np.median(h) > np.median(n)) == (mh > mn):
                agree += 1
        flag = "  <-- REPORTABLE" if agree >= MIN_CAMPAIGNS else ""
        print(f"{feat:<12}{mh:>13.2f}{mn:>12.2f}{agree:>10} of 4   "
              f"{'hits HIGHER' if mh > mn else 'hits LOWER':<12}{flag}")
        if agree >= MIN_CAMPAIGNS:
            keep.append((feat, mh, mn))

    print(f"\n(bar set before the numbers: a feature must agree on >= {MIN_CAMPAIGNS} of 4 "
          f"campaigns in the same direction)")
    if not keep:
        print("\n\U0001F534 NOTHING clears the bar. On these 82 variants no single structural")
        print("   property separates measured switches from non-switches consistently.")
    else:
        print(f"\n{len(keep)} feature(s) clear it: " + ", ".join(k[0] for k in keep))
        print("Per-campaign medians for those:")
        for feat, _, _ in keep:
            print(f"  {feat}")
            for name in cr.CAMPAIGNS:
                h = [r[feat] for r in rows if r["campaign"] == name and r["hit"]]
                n = [r[feat] for r in rows if r["campaign"] == name and not r["hit"]]
                if h and n:
                    print(f"     {name:<24}hits {np.median(h):>7.2f} (n={len(h)})   "
                          f"non {np.median(n):>7.2f} (n={len(n)})")
    json.dump(rows, open(os.path.join(HERE, "switch_features.json"), "w"), indent=1)
    print(f"\nwrote switch_features.json ({len(rows)} rows)")


if __name__ == "__main__":
    main()
