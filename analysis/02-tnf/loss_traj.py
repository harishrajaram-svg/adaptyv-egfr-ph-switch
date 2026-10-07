#!/usr/bin/env python3
"""Block-average a Mosaic run log's loss terms, to read the curve shape.

Mosaic re-predicts every step with a stochastic Boltz-2 pass, so single-step
values bounce by more than the whole optimisation moves. Block means are the
only honest way to ask whether a term is trending.

Usage:  loss_traj.py LOG [LOG ...] [--key .0.boltz2.0.bt_iptm] [--block 25]
        loss_traj.py --selftest
"""
import argparse
import re
import statistics
import sys

STEP = re.compile(r"^(\d+) loss:\s*(-?[0-9.]+)")


def parse(path, key):
    """Return (soft, sharp) lists of the term's value per step.

    The step counter restarts at 0 when the soft phase hands off to the sharp
    phase, which is the only marker of the boundary in the log.
    """
    phases, cur, last = [[]], None, -1
    cur = phases[0]
    pat = re.compile(re.escape(key) + r":\s*(-?[0-9.]+)")
    for line in open(path):
        m = STEP.match(line)
        if not m:
            continue
        step = int(m.group(1))
        if step < last:
            cur = []
            phases.append(cur)
        last = step
        v = pat.search(line)
        cur.append(float(v.group(1)) if v else None)
    soft = phases[0]
    sharp = phases[1] if len(phases) > 1 else []
    return soft, sharp


def blocks(vals, size):
    out = []
    for i in range(0, len(vals), size):
        x = [v for v in vals[i:i + size] if v is not None]
        out.append(statistics.mean(x) if x else None)
    return out


def selftest():
    """A restart in the step counter must split the phases, not concatenate."""
    import tempfile, os
    text = "".join(
        f"{s} loss: 1.0 .0.boltz2.0.bt_iptm: {v:.2f}\n"
        for s, v in list(enumerate([0.1] * 3)) + list(enumerate([0.9] * 2))
    )
    fd, p = tempfile.mkstemp(suffix=".log")
    os.write(fd, text.encode())
    os.close(fd)
    try:
        soft, sharp = parse(p, ".0.boltz2.0.bt_iptm")
        assert len(soft) == 3, soft
        assert len(sharp) == 2, sharp
        assert blocks(soft, 3) == [0.1]
        assert blocks(sharp, 2) == [0.9]
        # a missing key must not silently become 0.0
        soft2, _ = parse(p, ".0.boltz2.9.absent")
        assert blocks(soft2, 3) == [None], blocks(soft2, 3)
    finally:
        os.unlink(p)
    print("loss_traj.py --selftest PASS")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("logs", nargs="*")
    ap.add_argument("--key", default=".0.boltz2.0.bt_iptm")
    ap.add_argument("--block", type=int, default=25)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if not a.logs:
        ap.error("give at least one log, or --selftest")
    print(f"{a.key}  block={a.block}")
    for path in a.logs:
        soft, sharp = parse(path, a.key)
        f = lambda bs: " ".join("  -  " if b is None else f"{b:.3f}" for b in bs)
        print(f"  {path}")
        print(f"    soft  ({len(soft):3d}) {f(blocks(soft, a.block))}")
        if sharp:
            print(f"    sharp ({len(sharp):3d}) {f(blocks(sharp, a.block))}")


if __name__ == "__main__":
    sys.exit(main())
