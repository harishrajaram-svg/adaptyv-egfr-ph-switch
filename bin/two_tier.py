#!/usr/bin/env python3
"""Two-tier scoring: screen the whole pool at 1 seed, confirm the shortlist at 5.

Why two tiers. Max-over-seeds is a BIASED estimator -- it drifts upward with more
seeds -- so a 1-seed score and a 5-seed score are not comparable. Seeds buy you
less noise in the ORDERING, and that only changes decisions near the cut. Designs
at 0.89 and 0.00 order fine at one seed; designs at 0.62 and 0.61 do not order at
all. So spend the seeds where they change a decision.

Measured seed spread on barnase/barstar: 0.0099. Re-measure on a MARGINAL design.

Standing rule: scoring never eats the optimization rounds. If a run is behind,
cut seeds or sampling breadth before you cut predict-then-redesign cycles --
those are the cheapest hit-rate gain in the literature.

Usage:
  two_tier.py <designs_dir> --name <run> [--top 20] [--arm full|fast]
"""
import argparse, subprocess, shutil, sys, os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PY = ROOT / ".venv" / "bin" / "python"
ARMS = {"full": ("biohub/ESMFold2", "1afea82e432079d9af2ebd71d1e4c339ecca2ff0"),
        "fast": ("biohub/ESMFold2-Fast", "main")}


def fold(in_path: Path, out_dir: Path, run_name: str, seeds: str, arm: str, timeout: int):
    repo, rev = ARMS[arm]
    env = {**os.environ,
           "ESMFOLD2_HF_REPO": repo, "ESMFOLD2_HF_REVISION": rev,
           "MODAL_GPU": os.environ.get("MODAL_GPU", "L40S"),
           "MODAL_TIMEOUT": str(timeout),
           "PATH": f"{Path.home()}/.local/bin:" + os.environ.get("PATH", "")}
    cmd = ["modal", "run", "modal_esmfold2.py",
           "--input-faa", str(in_path.resolve()),
           "--seed", seeds,
           "--num-loops", "10", "--num-sampling-steps", "68",
           "--num-diffusion-samples", "1",
           "--out-dir", str(out_dir.resolve()), "--run-name", run_name]
    r = subprocess.run(cmd, cwd=ROOT / "biomodals", env=env)
    if r.returncode != 0:
        sys.exit(f"fold failed ({arm}, seeds={seeds}) -- exit {r.returncode}")


def rank(dir_: Path):
    """Return [(design, ipsae_min)] best-first, via bin/ipsae_min.py."""
    out = subprocess.run([str(PY), str(ROOT / "bin" / "ipsae_min.py"), "--dir", str(dir_)],
                         capture_output=True, text=True).stdout
    rows, seen = [], False
    for line in out.splitlines():
        if "MAX over seeds" in line:
            seen = True; continue
        if seen and "ipSAE_min=" in line:
            name = line.strip().split(":")[0]
            val = float(line.split("ipSAE_min=")[1].split()[0])
            rows.append((name, val))
    return sorted(rows, key=lambda r: -r[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("designs_dir")
    ap.add_argument("--name", required=True)
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--arm", default="full", choices=list(ARMS))
    ap.add_argument("--screen-seed", default="42")
    ap.add_argument("--final-seeds", default="1,2,3,4,5")
    a = ap.parse_args()

    src = Path(a.designs_dir)
    n_in = len(list(src.glob("*.faa")))
    base = ROOT / "runs" / a.name

    print(f"=== TIER 1: SCREEN — {n_in} designs, seed {a.screen_seed}, arm {a.arm} ===")
    fold(src, base / "screen", a.arm, a.screen_seed, a.arm, 60)
    screened = rank(base / "screen" / a.arm)
    if not screened:
        sys.exit("screen produced no scores")
    for n, v in screened[:a.top]:
        print(f"   {v:.4f}  {n}")

    keep = [n for n, _ in screened[:a.top]]
    print(f"\n=== TIER 2: CONFIRM — top {len(keep)}, seeds {a.final_seeds} ===")
    short = base / "shortlist"
    short.mkdir(parents=True, exist_ok=True)
    for n in keep:
        f = src / f"{n}.faa"
        if f.exists():
            shutil.copy(f, short / f.name)
    fold(short, base / "final", a.arm, a.final_seeds, a.arm, 60)
    final = rank(base / "final" / a.arm)

    sm = dict(screened)
    print(f"\n{'rank':<5}{'design':<34}{'5-seed':>9}{'1-seed':>9}{'Δ':>8}{'moved':>8}")
    order1 = {n: i for i, (n, _) in enumerate(screened)}
    for i, (n, v) in enumerate(final):
        d = v - sm.get(n, float('nan'))
        moved = order1.get(n, -1) - i
        print(f"{i+1:<5}{n:<34}{v:>9.4f}{sm.get(n, float('nan')):>9.4f}{d:>+8.4f}{moved:>+8d}")
    print("\nΔ is expected POSITIVE — max-over-seeds drifts up. Compare ranks, not raw scores,")
    print("across tiers. 'moved' = places gained vs the 1-seed order; large moves mean the")
    print("screen was noisy at this score range.")


if __name__ == "__main__":
    main()
