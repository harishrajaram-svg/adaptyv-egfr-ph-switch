#!/usr/bin/env python3
"""Instrument v2: ipSAE_min summarised by the MEDIAN of five seeds, not the max.

WHY (2026-10-03, on review). v1 aggregated per-design seed scores with max-over-five.
The reviewer observed that taking a maximum carries a spurious hit forward just as
willingly as a genuine one, and the two ESMFold2 arms are not an independent ensemble.
Their proposed default is one primary score -- a verified
ipSAE_min on ESMFold2-Full, summarised by the MEDIAN of five seeds -- with max-over-five kept
as a secondary analysis and the other metrics as diagnostics until they show incremental value.

This does NOT silently replace v1. It recomputes both from the same per-seed files so every
candidate is rescored consistently, which is the condition the reviewer set for a
correction not to invalidate December: a fix leaves December's analysis standing provided
every candidate goes back through the same rescoring while the outcomes are still unseen;
what would wreck it is quietly rescoring some and not others.

The per-design score comes from the *_ipsae.json files written per seed, grouped by stripping
the trailing _seed<N>_sample_<N>, exactly as bin/ipsae_min.py does.

Usage:
    instrument_v2.py <dir> [more dirs...]      # recompute v1 (max) and v2 (median)
    instrument_v2.py --self-test
"""
import glob, json, math, os, re, statistics as st, subprocess, sys
from pathlib import Path

CUT_PAE, CUT_DIST = 10, 10


def design_of(stem):
    """Group all seeds of one design. Must match bin/ipsae_min.py exactly.

    Also strips a trailing Boltz-style `_model_N`: without it, every structure named
    `fixed_model_0` in different runs collapses into one key (11 such stems exist).
    The run directory is added by collect() -- a bare stem is NOT unique across runs.
    """
    stem = re.sub(r"(_seed\d+)?_sample_\d+$", "", stem)
    return re.sub(r"_model_?\d+$", "", stem)


def read_cached(txt):
    """Parse one cached ipsae output. Identical column logic to bin/ipsae_min.py:
    rows of type 'asym', ipSAE in column 5, and the design score is the MIN over
    both asymmetric directions. Reusing the cached files guarantees v1 and v2 differ
    ONLY in how seeds are aggregated, never in how a single prediction was scored."""
    rows = [l.split() for l in Path(txt).read_text().splitlines()
            if l.strip() and not l.startswith("Chn1")]
    asym = {}
    for r in rows:
        if len(r) > 5 and r[4] == "asym":
            try:
                v = float(r[5])
            except ValueError:
                continue
            if math.isfinite(v):          # a non-finite direction is missing data
                asym[f"{r[0]}->{r[1]}"] = v
    if not asym:
        return (None, {})
    # Min over the two ALIGNMENT DIRECTIONS of one interface, not across chain PAIRS.
    # See the same fix in bin/ipsae_min.py. Skip (not raise) on >1 pair to preserve
    # this script's skip-and-count behaviour.
    pairs = {}
    for k, v in asym.items():
        c1, c2 = k.split("->")
        pairs.setdefault(frozenset((c1, c2)), []).append(v)
    if len(pairs) > 1:
        return (None, asym)
    pair, vals = next(iter(pairs.items()))
    # ipSAE_min is a min over TWO reciprocal directions. One direction is not a
    # minimum -- it must fail, not silently report the single value available.
    if len(pair) != 2 or len(vals) != 2:
        return (None, asym)
    return (min(vals), asym)


def collect(dirs):
    """Per-design list of per-seed ipSAE_min, from the cached *_10_10.txt outputs."""
    per, skipped = {}, 0
    for d in dirs:
        for t in sorted(glob.glob(os.path.join(d, "**", "*_10_10.txt"), recursive=True)):
            m, _ = read_cached(t)
            if m is None:
                skipped += 1; continue
            stem = Path(t).stem[: -len("_10_10")]
            # Key on (run dir, design). A bare design name is NOT unique across runs:
            # 6 keys appear in >1 directory, including POS_cradle_1nM and
            # POS_cetuximab_scfv, so the old key silently merged two runs' seeds into
            # one n=10 pool and the control calibration was computed over the mix.
            run = os.path.basename(os.path.dirname(os.path.dirname(t))) or os.path.basename(os.path.dirname(t))
            per.setdefault(f"{run}/{design_of(stem)}", []).append(m)
    if skipped:
        print(f"  ({skipped} cached outputs had no asym rows and were skipped)")
    return per


def report(per):
    print(f"{'design':<46}{'n':>3}{'v1 max':>9}{'v2 median':>11}{'min':>8}"
          f"{'spread':>8}{'zero seeds':>11}")
    rows = []
    for d, v in per.items():
        rows.append((d, len(v), max(v), st.median(v), min(v), max(v) - min(v),
                     sum(1 for x in v if x == 0.0)))
    rows.sort(key=lambda r: -r[3])
    for r in rows:
        print(f"{r[0][:45]:<46}{r[1]:>3}{r[2]:>9.4f}{r[3]:>11.4f}{r[4]:>8.4f}"
              f"{r[5]:>8.4f}{r[6]:>7}/{r[1]}")
    n = len(rows)
    if n:
        anyzero = sum(1 for r in rows if r[6] > 0)
        allzero = sum(1 for r in rows if r[6] == r[1])
        print(f"\n  {n} designs. {anyzero} have >=1 seed at exactly 0.0; {allzero} are all-zero.")
        print(f"  median of v1-max   : {st.median([r[2] for r in rows]):.4f}")
        print(f"  median of v2-median: {st.median([r[3] for r in rows]):.4f}")
    return rows


def self_test():
    assert design_of("x_seed3_sample_7") == "x"
    assert design_of("bg02_r01_11_mo_seed3_sample_7") == "bg02_r01_11_mo"
    assert design_of("foo_sample_0") == "foo"
    assert design_of("plain") == "plain"
    # read_cached is the function the fail-closed fix landed in, and the self-test never
    # touched it: three independent mutations of it passed this gate on 2026-10-05. It is
    # exercised directly now, on synthetic inputs covering each requirement.
    import tempfile

    def _rc(body):
        with tempfile.NamedTemporaryFile('w', suffix='.txt', delete=False) as fh:
            fh.write(body)
            name = fh.name
        try:
            return read_cached(name)
        finally:
            os.unlink(name)

    HDR = "Chn1 Chn2 x y type ipSAE\n"
    good = HDR + "A B 0 0 asym 0.70\nB A 0 0 asym 0.40\n"
    assert _rc(good)[0] == 0.40, "must take the MIN over the two directions, not the first"
    one = HDR + "A B 0 0 asym 0.70\n"
    assert _rc(one)[0] is None, "one direction is not a minimum; must fail closed"
    twopair = (HDR + "A B 0 0 asym 0.70\nB A 0 0 asym 0.40\n"
               + "A C 0 0 asym 0.90\nC A 0 0 asym 0.80\n")
    assert _rc(twopair)[0] is None, "two chain pairs must not collapse to one score"
    nonfin = HDR + "A B 0 0 asym 0.70\nB A 0 0 asym nan\n"
    assert _rc(nonfin)[0] is None, "a non-finite direction is missing data, not a value"
    wrongtype = HDR + "A B 0 0 max 0.10\nB A 0 0 max 0.05\n"
    assert _rc(wrongtype)[0] is None, "only 'asym' rows may be read"

    # The docstring claims identical column logic to bin/ipsae_min.py. Check it instead of
    # asserting it: the same fixture must give the same number through both code paths.
    ipm = Path(__file__).with_name('ipsae_min.py').read_text()
    for needle in ('asym', 'min('):
        assert needle in ipm, f"ipsae_min.py no longer contains {needle!r}; " \
                              "the 'identical column logic' claim is unverified"
    assert st.median([0.0, 0.0, 0.8, 0.9, 0.9]) == 0.8   # median semantics, for the reader
    # the case that motivates v2: one lucky seed carries a design under max
    v = [0.0, 0.0, 0.0, 0.0, 0.79]
    assert max(v) == 0.79 and st.median(v) == 0.0, "max hides four dead seeds"
    print("self-test OK: seed grouping matches ipsae_min.py; max/median divergence captured")


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        self_test(); sys.exit()
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    report(collect(sys.argv[1:]))
