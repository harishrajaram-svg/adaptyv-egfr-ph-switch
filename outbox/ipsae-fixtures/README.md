# ipSAE_min — fixture package for independent review

Built for PK, 2026-10-04, in answer to: *"please send the actual file and reproducible
fixtures ... include the code commit, model/version, raw PAE and structure files,
chain/residue mapping, cutoffs, and expected directional scores."*

## Run it

    python3 run_fixtures.py            # score every case, both directions, then the min
    python3 run_fixtures.py --check    # compare to expected.json, exit 1 on mismatch

No arguments, no network. Needs `gemmi` and whatever `ipsae/ipsae.py` imports (numpy).

## Provenance

| | |
|---|---|
| our wrapper | `bin/ipsae_min.py` (copied here as `ipsae_min.py`) |
| repo commit | `76c5b0c0920fde7be2f4b8324ca9be1db8b0f304`, 2026-10-04 00:42 EDT, **with `bin/ipsae_min.py` modified in the working tree** — the copy here is the working-tree version, not the committed one |
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

## The eight cases

| case | construct | chains | ipSAE_min | why it is here |
|---|---|---|---|---|
| 02 | full ECD | A=621 target, B=127 VHH | **0.000000** | a real 5 nM nanobody scoring zero |
| 03 | d3 crop | A=170 target, B=127 VHH | **0.000000** | same molecule, cropped target |
| 04 | full ECD | A=621, B=241 scFv | 0.639298 | real positive |
| 05 | d3 crop | A=170, B=241 scFv | **0.400501** | **same molecule, different construct** |
| 06 | d3 crop | **A=131 binder**, B=170 target | 0.000000 | chain order swapped vs the design pool |
| 07 | full ECD | A=621, B=134 | 0.140240 | the EGF-derived control, activity unknown |
| 08 | d3 crop | A=150 binder, B=170 target | 0.576505 | rank 1 |

Barnase/barstar is not in the package — the cached pair is not on disk in a form this runner
can reproduce. Treat the 0.8887 regression number as unverified until it is rebuilt.

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
`ipTM_af` is identically **0.000** in all 4,192 cached outputs, and `pDockQ` / `pDockQ2` are
constant at **0.0183 / 0.0073**. Cause: our ESMFold2 sidecar writes pLDDT on a 0–1 scale into
the AF2 code path, which rescales only on the Boltz branch (`ipsae.py:472`). Nothing we rank
on uses them.

## Still open on our side

- barnase/barstar regression rebuilt from source files
- an empty-interface case where the binder is deliberately placed far from the target
- a synthetic 3-chain complex to exercise the >1-pair refusal path
