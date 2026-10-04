# ipSAE_min — fixture package for independent review

Built for PK, 2026-10-04, in answer to: *"please send the actual file and reproducible
fixtures ... include the code commit, model/version, raw PAE and structure files,
chain/residue mapping, cutoffs, and expected directional scores."*

## Run it

    python3 run_fixtures.py            # score every case, both directions, then the min
    python3 run_fixtures.py --check    # compare to expected.json, exit 1 on mismatch
    python3 run_fixtures.py --freeze   # re-freeze expected.json (deliberate only)

**Self-contained.** No arguments, no network, nothing outside this directory. Needs `gemmi` and
`numpy`. The pinned Dunbrack reference is vendored at `vendor/ipsae/ipsae.py` under its MIT
licence (commit `6174cf9e71cb1bd660cc805856a18c4871a6dec3`), and the runner prefers a repo
checkout at `../../ipsae/` when one exists so the two copies can be diffed.

An earlier version of this bundle was **not** runnable on your machine: it scored through
`../../ipsae/ipsae.py`, and `ipsae/` is gitignored in our repo, so standalone every case errored
out. Vendoring fixes that.

`--check` compares against `expected.json`, which is now frozen: a plain run no longer
overwrites it. It used to, which means a regression would silently have become the new
expectation — the same error as choosing a threshold after seeing the data.

Expected output: **10/10 cases reproduce**, exit 0.

## Provenance

| | |
|---|---|
| our wrapper | `bin/ipsae_min.py` (copied here as `ipsae_min.py`) |
| repo commit | `a5bd3d6`, 2026-10-04. `bin/ipsae_min.py` is **committed and the tree is clean**, so the copy here is byte-identical to the committed one. (An earlier version of this table pinned `76c5b0c` and warned the scorer was dirty in the working tree; that caveat no longer applies.) |
| public repo | <https://github.com/harishrajaram-svg/adaptyv-egfr-ph-switch> — this bundle, the methods document and the submission are all there |
| reference implementation | `ipsae/ipsae.py`, Roland Dunbrack, version 4, biorxiv 2025.02.10.637595v2 |
| reference commit | `6174cf9e71cb1bd660cc805856a18c4871a6dec3`, 2026-01-03 (shallow clone) |
| structure model | ESMFold2 via Modal, Anthropic's published protocol parameters (loops=10, steps=68) |
| cutoffs | PAE 10, distance 10 Å — passed explicitly on every call |
| aggregation under review | min over the two **alignment directions** of one interface |

## What we compute, and the distinction you flagged

The reference emits one row per direction (`A->B`, `B->A`) plus its own `max` row. We take
the **minimum of the two asymmetric directions**, never the reference's `max`, and never a
minimum across chain *pairs* — on a 3-chain complex that would let an intra-binder interface
become the score. `bin/ipsae_min.py` raises rather than guessing when >1 pair is present.

**This is a different question from max-versus-median across seeds.** That one is seed
aggregation and lives in `bin/instrument_v2.py`; your max→median correction applies there, and
we hold 5 seeds with the median as the provisional primary summary, retaining per-seed values.

Column mapping (reference output, 0-indexed): `5` = ipSAE (what we use), **`13` = n0res, the
count of residue pairs passing the interface filter**, `14` = n0chn, the total chain-pair
residue count. An earlier version of this script printed column 14 and labelled it "interface
residues", which made a zero-interface case look like a full interface scoring zero. Fixed.

## The ten cases

| case | construct | chains | ipSAE_min | why it is here |
|---|---|---|---|---|
| 01 | barnase/barstar | A=110, B=89 | **0.888672** | the positive-control regression you named first. A→B 0.888672, B→A 0.939795 |
| 02 | full ECD | A=621 target, B=127 VHH | **0.000000** | a real 5 nM nanobody scoring zero |
| 03 | d3 crop | A=170 target, B=127 VHH | **0.000000** | same molecule, cropped target |
| 04 | full ECD | A=621, B=241 scFv | 0.639298 | real positive |
| 05 | d3 crop | A=170, B=241 scFv | **0.400501** | **same molecule, different construct** |
| 06 | d3 crop | **A=131 binder**, B=170 target | 0.000000 | chain order swapped vs the design pool |
| 07 | full ECD | A=621, B=134 | 0.140240 | the EGF-derived control, activity unknown |
| 08 | d3 crop | A=150 binder, B=170 target | 0.576505 | a submitted design (`rimA01_r15`, now rank 6 of 10) |
| 09 | full ECD | A=159 shuffled, B=621 target | **0.000000** | **empty interface**: 0 of 780 residue pairs pass the filter, in BOTH directions |
| 10 | d3 crop | A=204 target, B+C = Fv heavy/light | **refuses** | 3 inter-chain pairs. `ipSAE_min` raises rather than letting the intra-Fv B:C interface (0.8659) become the binder score |

**Correction, 2026-10-04.** An earlier version of this README told you barnase/barstar was
"not on disk in a form this runner can reproduce" and to treat 0.8887 as unverified. That was
wrong. The cached pair was on disk the whole time, in exactly the `*.cif` / `*_ipsae.json`
pattern this runner globs, and it reproduces to **0.888672** — the number you quoted back to
me. It is now **case 01**, and it is the bundle's only true-positive regression. There is also
a five-seed set at `runs/seedtest/pos5/` that nobody had found: 0.8800–0.8899, median 0.8878,
spread 0.0099, which is the seed-level evidence for keeping the median as the primary summary.

## The VHH zeros, traced end to end

You asked that a failed run must not silently become a valid score of zero. It does not.
Traced through parse → PAE → reference output → interface mask → aggregation:

```
case 02  POS_nano2_5nM (real 5 nM VHH), full ECD
  1 parse       chains A=621, B=127          OK
  2 PAE         748x748, range 0.25-31.67    OK, present and sane
  3 reference   3 rows, 2 asym               OK, no exception, no missing file
  4 interface   n0res = 0 in BOTH directions <-- the cause
  5 aggregate   min(0.000000, 0.000000) = 0
```

**Zero residue pairs pass the PAE<10 / distance<10 interface filter.** The model predicts no
confident interface at all. That is a genuine model output, not a crash: a failure would show
a missing PAE, no output file, or an exception, and none occurs.

Case 03 is the informative companion — the same VHH against the cropped target gives
`A->B` 0.0147 with **14** interface residues and `B->A` 0.0000 with **0**. The directional
asymmetry is real, and our min-over-directions rule takes the zero.

So: the pipeline does not mis-handle nanobodies, ESMFold2 does not predict their interfaces
here. Per your instruction we now mark VHH-format designs **inadequately assessed by this
pipeline** rather than rejecting them on score.

## Construct-dependent calibration, in one comparison

Cases 04 and 05 are **the same cetuximab scFv sequence, the same reference, the same cutoffs**,
differing only in target construct:

    full ECD (621 aa target)   ipSAE_min 0.639
    domain III crop (170 aa)   ipSAE_min 0.401

A 0.24 shift from the construct alone. Any threshold must come from controls folded against
the same construct as the design it judges. We previously applied one ECD-derived bar to both.

## Known-invalid columns in the reference output

Not used by us, flagged so they are not reached for as an orthogonal check:
`ipTM_af` is **0.000** in 18,006 of the 18,039 interface rows across our 6,011 cached outputs,
and `pDockQ` / `pDockQ2` are constant at **0.0183 / 0.0073**. The 33 exceptions are Boltz
outputs, which take the rescale branch. An earlier version of this file said "identically 0.000
in all 4,192 cached outputs", which was wrong on the count and on the universality. Cause: our ESMFold2 sidecar writes pLDDT on a 0–1 scale into
the AF2 code path, which rescales only on the Boltz branch (`ipsae.py:472`). Nothing we rank
on uses them.

## Your named checks, and where each one is

You asked for: *"an asymmetric complex, chain-order swaps, full-ECD versus cropped inputs, and
empty-interface cases"*, plus the barnase/barstar regression, plus a trace of the VHH zeros.

| your check | case | what it shows |
|---|---|---|
| barnase/barstar regression | **01** | 0.888672, matching the number you quoted |
| VHH zeros traced end to end | **02** | 0 interface residues of 748, both directions — absent, not a failed run |
| asymmetric complex | **03** | A→B 0.014687 vs B→A 0.000000 — a 0.0147 directional split on one interface, which is why the min is taken over directions |
| full ECD vs cropped target | **04 / 05** | same scFv, same cutoffs: 0.639 on the 621 aa ECD, 0.401 on the 170 aa crop. A 0.24 shift from the construct alone |
| chain-order swap | **06** | binder in chain A rather than B; 29 of our 69 run directories have the target first |
| empty interface | **09** | a composition-matched shuffled null: 0 interface residues of 780 in **both** directions |
| >1-pair refusal path | **10** | a 3-chain Fv complex. `ipSAE_min` **refuses** rather than letting an intra-Fv interface become the binder score. `--check` asserts the refusal, so "the right answer is an error" is a testable expectation |

## Still open on our side

- **Your independent review of the scorer itself.** This bundle demonstrates behaviour; it does
  not substitute for a second pair of eyes on the 80 lines. That remains the ask.
- A second structure model. Only ESMFold2-Full has ever scored a design here: ESMFold2-Fast,
  which you specified as a robustness check, has four controls and no designs. Chai-1 and
  Protenix v2 were never built.
- No case yet covers an indexing error in the residue mapping, which you raised separately from
  chain order.
