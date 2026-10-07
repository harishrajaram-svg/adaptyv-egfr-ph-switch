#!/usr/bin/env python3
"""Split a Mosaic design wave into jobs that fit Modal's function ceiling, and print the plan.

WHY THIS EXISTS. The wave as specified does not fit in one job. At the measured rate
(30.6 s/step for the three-leg human-trimer + mouse-dimer configuration, s38) 16 trajectories
x 200 steps is ~33 hours against Modal's typical 24 h function timeout. It has to be split.

WHY IT ONLY PRINTS. Launching is behind --launch, off by default, because a wave costs real
money and needs an explicit decision. `plan` is safe to run and re-run.

HOW IT SPLITS. By LENGTH GROUP, keeping every seed for a length in the same job. The CLI can
only express "these lengths x seeds seed0..seed0+n", so length groups are the natural unit --
and it means each job is independently interpretable: a finished job is a complete answer for
its lengths rather than a fragment of every length.

RATE, AND WHY IT CARRIES A FACTOR. The same configuration measured 25.9 s/step on one run and
30.6 on another -- identical tokens, legs and length, 18% apart -- and the ESMFold estimator
missed by a similar margin in the other direction. The default rate is the slower observation
with a 1.15 factor. Pass --s-per-step to override from a fresh smoke; do NOT lower the default
to a point estimate.

usage:
  python3 bin/plan_wave.py --lengths 68,72,76,80,84,88 --n-seeds 2 --steps 200
  python3 bin/plan_wave.py ... --ceiling-hours 12 --launch
  python3 bin/plan_wave.py --selftest
"""
import math
import shlex
import subprocess
import sys

DEFAULT_S_PER_STEP = 30.6 * 1.15      # see the note above; slower observation + factor
COMPILE_S = 703.0                      # measured, p2dry03 -- per job, not per trajectory
REPRED_S = 300.0                       # step-4 re-prediction, per trajectory
MODAL_CEILING_H = 24.0                 # Modal's typical function timeout


def job_minutes(n_traj, steps, s_per_step=DEFAULT_S_PER_STEP):
    """Wall clock for one job: one compile, then every trajectory's steps + re-prediction."""
    return (COMPILE_S + n_traj * (steps * s_per_step + REPRED_S)) / 60.0


def split(lengths, n_seeds, steps, ceiling_h, s_per_step=DEFAULT_S_PER_STEP):
    """Group lengths so each job fits `ceiling_h`. Returns [[lengths], ...].

    Fails closed: if a SINGLE length with all its seeds cannot fit, that is not a splitting
    problem and the caller must cut steps or seeds instead -- so it raises rather than
    emitting a job it knows will be truncated.
    """
    ceiling_min = ceiling_h * 60.0
    one = job_minutes(n_seeds, steps, s_per_step)
    if one > ceiling_min:
        raise SystemExit(
            f"REFUSING TO PLAN: one length x {n_seeds} seeds x {steps} steps is "
            f"{one:.0f} min, over the {ceiling_min:.0f} min ceiling. Splitting cannot fix "
            f"this -- reduce --steps or --n-seeds, or raise --ceiling-hours."
        )
    per_job = max(1, int(ceiling_min // one))          # lengths that fit in one job
    return [lengths[i:i + per_job] for i in range(0, len(lengths), per_job)]


def commands(groups, n_seeds, steps, run_prefix, passthrough, soft_frac=0.75):
    """One exact command per job. Soft:sharp keeps the wave's 3:1 ratio."""
    soft = int(round(steps * soft_frac))
    out = []
    for k, grp in enumerate(groups, start=1):
        out.append(
            f"GPU=L40S TIMEOUT={math.ceil(job_minutes(len(grp) * n_seeds, steps)) + 15} "
            f"modal run --detach modal_mosaic.py --step 4 "
            f"--run-name {run_prefix}-j{k} "
            f"--lengths {','.join(str(v) for v in grp)} "
            f"--n-seeds {n_seeds} --steps-soft {soft} --steps-sharp {steps - soft} "
            + passthrough
        )
    return out


def selftest():
    # a job that fits is not split further than it needs to be
    g = split([76], 2, 200, 24.0)
    assert g == [[76]], g
    # six lengths at 200 steps must split, and every job must fit the ceiling
    g = split([68, 72, 76, 80, 84, 88], 2, 200, 12.0)
    assert len(g) > 1, g
    assert sum(len(x) for x in g) == 6 and [v for x in g for v in x] == [68, 72, 76, 80, 84, 88], g
    for grp in g:
        assert job_minutes(len(grp) * 2, 200) <= 12 * 60 + 1e-6, (grp, job_minutes(len(grp) * 2, 200))
    # nothing is dropped or duplicated for any input size
    for n in range(1, 13):
        gg = split(list(range(n)), 2, 200, 12.0)
        assert [v for x in gg for v in x] == list(range(n)), (n, gg)
    # MUTATION TEST: a ceiling that cannot fit even one length must RAISE, not emit a job
    # that will be truncated. Silently emitting it is how a wave gets billed for nothing.
    try:
        split([76], 2, 200, 0.5)
        raise AssertionError("an impossible ceiling must raise, not plan a truncated job")
    except SystemExit:
        pass
    # the wave, as the numbers actually stand
    one_job = job_minutes(16, 200)
    assert one_job / 60 > MODAL_CEILING_H, (
        f"one-job wave is {one_job/60:.1f} h -- if this no longer exceeds {MODAL_CEILING_H} h, "
        f"the splitting rationale has changed and the comments need revisiting")
    # every emitted command carries a timeout ABOVE its own estimate
    cmds = commands(g, 2, 200, "t", "--target x")
    for grp, c in zip(g, cmds):
        t = int([p for p in c.split() if p.startswith("TIMEOUT=")][0].split("=")[1])
        assert t >= job_minutes(len(grp) * 2, 200), (t, grp)
    assert all("--run-name t-j" in c for c in cmds)
    assert len({c.split("--run-name ")[1].split()[0] for c in cmds}) == len(cmds), \
        "run names must be unique or jobs overwrite each other on the Volume"
    print(f"plan_wave     one-job wave = {one_job/60:.1f} h, over the {MODAL_CEILING_H} h "
          f"ceiling, so splitting is required")
    print(f"              6 lengths/2 seeds/200 steps -> {len(g)} jobs, each within 12 h; "
          f"no length dropped or duplicated for n=1..12")
    print("              impossible ceiling raises; every TIMEOUT exceeds its own estimate; "
          "run names unique")
    print("selftest OK")


def main():
    a = sys.argv[1:]
    if "--selftest" in a:
        return selftest()
    g = lambda f, d: a[a.index(f) + 1] if f in a else d
    lengths = [int(v) for v in g("--lengths", "76").split(",")]
    n_seeds = int(g("--n-seeds", "2"))
    steps = int(g("--steps", "200"))
    ceiling = float(g("--ceiling-hours", "12"))
    sps = float(g("--s-per-step", str(DEFAULT_S_PER_STEP)))
    prefix = g("--run-prefix", "p2wave")
    # everything after --passthrough is handed to every job verbatim
    pt = " ".join(a[a.index("--passthrough") + 1:]) if "--passthrough" in a else ""

    groups = split(lengths, n_seeds, steps, ceiling, sps)
    total_traj = len(lengths) * n_seeds
    print(f"wave: {len(lengths)} lengths x {n_seeds} seeds = {total_traj} trajectories "
          f"x {steps} steps at {sps:.1f} s/step")
    print(f"one job would be {job_minutes(total_traj, steps, sps)/60:.1f} h "
          f"(Modal ceiling ~{MODAL_CEILING_H} h) -> {len(groups)} jobs\n")
    cmds = commands(groups, n_seeds, steps, prefix, pt)
    grand = 0.0
    for k, (grp, c) in enumerate(zip(groups, cmds), start=1):
        mins = job_minutes(len(grp) * n_seeds, steps, sps)
        grand += mins
        print(f"  job {k}: lengths {grp}  {len(grp)*n_seeds} traj  "
              f"~{mins:.0f} min ({mins/60:.1f} h)")
    print(f"\ntotal ~{grand/60:.1f} GPU-hours, ~${grand/60*1.95:.2f} on L40S")
    print("jobs are independent: each writes its own run-name on the Volume, and designs are "
          "persisted per-design,\nso a failed job costs only itself and can be re-run alone.\n")
    for c in cmds:
        print(c + "\n")
    if "--launch" not in a:
        print("DRY PLAN. Nothing launched. Add --launch to run these sequentially.")
        return
    for k, c in enumerate(cmds, start=1):
        print(f"[plan_wave] launching job {k}/{len(cmds)}")
        r = subprocess.run(c, shell=True, cwd="biomodals")
        if r.returncode != 0:
            raise SystemExit(f"job {k} exited {r.returncode} -- stopping; later jobs NOT "
                             f"launched. Earlier jobs' designs are on the Volume.")


if __name__ == "__main__":
    main()
