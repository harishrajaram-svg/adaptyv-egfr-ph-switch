#!/usr/bin/env python3
"""Score the free-footprint probe against the pinned probe, and name the verdict.

One command, because the reading that matters is easy to get wrong: the decision
turns on the TRAJECTORY SHAPE, not the final median. Reading only endpoints is what
made the 200-step test cost $15 for half an answer (METHODS s11).

Usage:  score_free_probe.py                      # once the four runs have landed
        score_free_probe.py --selftest
"""
import argparse
import csv
import math
import os
import re
import statistics
import sys

RUNS = "runs/mosaic-p2"
PINNED = {("76", "0"): "p2probe-a", ("76", "1"): "p2probe-b",
          ("84", "0"): "p2probe-c", ("84", "1"): "p2probe-d"}
FREE = {("76", "0"): "p2free-p", ("76", "1"): "p2free-q",
        ("84", "0"): "p2free-r", ("84", "1"): "p2free-s"}
LOGS = {"p2free-p": "/tmp/free_p.log", "p2free-q": "/tmp/free_q.log",
        "p2free-r": "/tmp/free_r.log", "p2free-s": "/tmp/free_s.log"}
KEY = ".0.boltz2.0.bt_iptm"
STEP = re.compile(r"^(\d+) loss:")


def row(run):
    f = os.path.join(RUNS, run, "designs.tsv")
    if not os.path.isfile(f):
        return None
    with open(f) as fh:
        rs = list(csv.DictReader(fh, delimiter="\t"))
    return rs[0] if rs else None


def soft(path):
    """The soft-phase values of KEY. The step counter restarts at the sharp handoff."""
    out, last = [], -1
    pat = re.compile(re.escape(KEY) + r":\s*(-?[0-9.]+)")
    if not os.path.isfile(path):
        return out
    for line in open(path):
        m = STEP.match(line)
        if not m:
            continue
        if int(m.group(1)) < last:
            break
        last = int(m.group(1))
        v = pat.search(line)
        out.append(float(v.group(1)) if v else None)
    return [v for v in out if v is not None]


def blocks(vals, size):
    return [statistics.mean(v) for i in range(0, len(vals), size)
            if (v := [x for x in vals[i:i + size] if x is not None])]


def slope(y):
    """OLS slope per step and its standard error. None when too short to fit."""
    n = len(y)
    if n < 10:
        return None, None
    mx, my = (n - 1) / 2, sum(y) / n
    sxx = sum((i - mx) ** 2 for i in range(n))
    b = sum((i - mx) * (v - my) for i, v in enumerate(y)) / sxx
    a = my - b * mx
    s2 = sum((v - (a + b * i)) ** 2 for i, v in enumerate(y)) / (n - 2)
    return b, math.sqrt(s2 / sxx)


def selftest():
    assert blocks([1.0, 3.0, 2.0, 4.0], 2) == [2.0, 3.0]
    assert blocks([], 2) == []
    assert blocks([1.0, None, 3.0], 3) == [2.0]          # None skipped, not zeroed
    assert slope([1.0] * 5) == (None, None)              # too short to fit
    b, se = slope([float(i) for i in range(20)])
    assert abs(b - 1.0) < 1e-9 and se < 1e-9, (b, se)    # a clean ramp recovers slope 1
    b, se = slope([5.0] * 20)
    assert abs(b) < 1e-12, b                             # flat recovers slope 0
    print("score_free_probe.py --selftest PASS")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return selftest()

    print("PAIRED ENDPOINTS  (same length, same seed, same weights; footprint is the variable)")
    print(f"  {'len/seed':>9} {'pinned 50':>10} {'free 100':>9} {'diff':>8}")
    pin, fre, missing = [], [], []
    for k in sorted(PINNED):
        rp, rf = row(PINNED[k]), row(FREE[k])
        if rf is None:
            missing.append(FREE[k])
            continue
        vp, vf = float(rp["iptm_repred"]), float(rf["iptm_repred"])
        pin.append(vp); fre.append(vf)
        print(f"  {k[0]+' / '+k[1]:>9} {vp:>10.4f} {vf:>9.4f} {vf-vp:>+8.4f}")
    if missing:
        print(f"\n  NOT LANDED YET: {', '.join(missing)}")
    if fre:
        print(f"\n  {'mean':>9} {statistics.mean(pin):>10.4f} {statistics.mean(fre):>9.4f} "
              f"{statistics.mean(fre)-statistics.mean(pin):>+8.4f}")
        print(f"  {'median':>9} {statistics.median(pin):>10.4f} {statistics.median(fre):>9.4f}")
        print(f"  {'max':>9} {max(pin):>10.4f} {max(fre):>9.4f}")

    print("\nTRAJECTORY SHAPE  <-- THE READING THAT DECIDES  (25-step block means)")
    slopes = []
    for k in sorted(FREE):
        y = soft(LOGS[FREE[k]])
        if not y:
            continue
        bs = blocks(y, 25)
        b, se = slope(y)
        if b is not None:
            slopes.append((b, se))
            rise, mde = b * len(y), 2 * se * len(y)
            verdict = "CLIMBS" if rise > mde else ("falls" if -rise > mde else "flat")
            print(f"  {FREE[k]:<10} n={len(y):>3}  " + " ".join(f"{v:.3f}" for v in bs)
                  + f"   rise {rise:+.3f} +/- {mde:.3f}  {verdict}")
        else:
            print(f"  {FREE[k]:<10} n={len(y):>3}  too short to fit")
    if len(slopes) == 4:
        tot = sum(b for b, _ in slopes) / 4
        n = len(soft(LOGS["p2free-p"]))
        mde = 2 * math.sqrt(sum(se ** 2 for _, se in slopes)) / 4 * n
        rise = tot * n
        print(f"\n  FOUR-RUN MEAN rise over the soft phase: {rise:+.4f} +/- {mde:.4f}")
        print("\nVERDICT")
        if rise > mde:
            print("  The interface term CLIMBS with a free footprint. The pinned epitope was")
            print("  the blocker. -> ROADMAP item 10: the wave is worth ~$49, weighted to free.")
        else:
            print("  FLAT, like all 12 pinned trajectories. The cause is not the footprint:")
            print("  it is gradient noise or term interaction (METHODS s6.3 causes 2 and 3).")
            print("  -> plan C, the methods-first submission. Nothing reachable in the")
            print("     remaining time clears 0.45.")
    else:
        print("\n  (verdict needs all four runs)")


if __name__ == "__main__":
    sys.exit(main())
