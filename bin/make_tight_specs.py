#!/usr/bin/env python3
"""Arms built on the one defensible switch of 2026-10-03, and on why it still failed.

WHAT WE LEARNED TONIGHT, and what each arm does about it.

g532mimic_short_14 is the first design all week with a real pH switch: H433 pKa
6.37 -> 7.27 (+0.90, 2.48x stronger at pH 6.5), from a PINNED Asp21 sitting 2.72 A
from the histidine ring. A knockout control (D21->Ala, nothing else changed) puts
the causal effect of that single carboxylate at +2.09 pKa units, and shows the
burial-only case does exactly what the week's negative predicted: -1.19, acid
weakens. So the negative was never wrong, only incomplete -- burial does suppress
the switch, and one carboxylate close enough overcomes the penalty and flips it.

Three facts then decide the arms:

  1. PLACEMENT, not identity, is the bottleneck. Pins hold 80/80, but the
     salt-bridge geometry (<4.0 A) occurred in only 2/160 design-histidine pairs.
     Cause: `binding: 343,406` let the interface spread over 24 target residues.
     -> ARM tight_h433 narrows the declared target site to H433 +/- 3 residues
        (positional 403..409 = T430 K431 Q432 H433 G434 Q435 F436, which is the
        RTKQHGQF pose-check motif minus its Arg -- an independent numbering check).

  2. H370 IS DEAD WEIGHT. 61% buried in apo, pKa 5.53, median dPKa -0.58 across
     80 designs, nothing positive. It consumed half of every mimic spec.
     -> dropped from every arm here, which doubles the useful sampling.

  3. H433 IS CAPPED AT 5.41x. It has no exposed histidine partner within 25 A, so
     it is a single-site target and that is its thermodynamic ceiling. The only
     place on EGFR where a >10x switch is geometrically reachable is a domain IV
     cluster: canonical H584 + H615, 6.9 A apart, combined ceiling 23.4x.
     -> ARM dom4_pair. Caveat recorded, not hidden: mature 560 sits in the tether
        region, so this collides with the conformational question the reviewer raised for
        domain II. High information either way.

The carboxylate arms vary only the pinned chemistry at the tight site, holding the
binder band and the target patch fixed, so dPKa is the only thing that moves.

CONSTRAINTS, all learned by breaking them:
  * every a..b range must come AFTER a pinned residue, or its index is re-sampled
    and no label can address it;
  * a loop window must not cover a pinned residue -- BoltzGen requires every
    SS-labelled residue to be designed, and a pinned letter is not (data.py:1953);
  * binder `binding_types` must be EXACTLY the pinned indices. A binding type on a
    designed residue is rejected (data.py:1949); the rim arms died of this;
  * letters must be UPPERCASE -- lowercase silently becomes UNK.
"""
import sys
from pathlib import Path

# positional indices within egfr_ecd_6aru_renum.pdb (contiguous 1..609, renum = canonical - 27)
TGT_H433 = "403,404,405,406,407,408,409"   # canonical 430..436, H433 at the centre
TGT_DOM4 = "555,556,557,558,559,586,587,588,589,590"  # canonical H584 (renum 557) + H615 (renum 588), each +/- 2

# (sequence, pinned indices, loop windows)  -- binder band 56..90, the only band that hit
BAND = (56, 90)   # the only binder band that produced a switch; 123-160 and 31-49 gave zero

BINDER = {
    "DE": ("20, D, 14, E, 20..54", "21,36", "16..20,22..26,31..35,37..41"),  # the winner's own spec
    "DD": ("20, D, 14, D, 20..54", "21,36", "16..20,22..26,31..35,37..41"),
    "EE": ("20, E, 14, E, 20..54", "21,36", "16..20,22..26,31..35,37..41"),
    "D1": ("20, D, 35..69",        "21",    "16..20,22..26"),
}

ANCHORS = {406, 557, 588,          # renum of canonical H433, H584, H615 -- verified HIS
            123}                   # H433's POSITIONAL index in the domain-III crop (406 - 284 + 1)

# The crop arms. Added 2026-10-04 after the Mosaic session's 5/5 epitope result turned out to
# be explained by its domain-III CROP rather than by its optimized contact objective: our own
# conditional rate (among designs that bound domain III at all) is 54/65 = 83% patch contact,
# against its 5/5, Fisher p = 1.000. The crop is doing the work, and a crop is free.
#
# We ran every arm against the FULL 609-residue ECD, so about half of every arm escaped to
# domain I. Cropping removes that escape route with no new tooling. This is the clean A/B:
# identical binder spec, identical declared patch, CROP IS THE ONLY CHANGE vs tight_h433
# (which gave 0/30 geometry and only 6/30 majority-domain-III).
#
# CROP RENUMBERS EVERYTHING. In egfr_d3_crop_renum.pdb (renum 284..453, 170 residues) the
# anchor H433 is positional 123, not 406, and the patch is 120..126. Verified residue by
# residue: positional 120..126 = THR LYS GLN HIS GLY GLN PHE. The Mosaic session independently
# reports anchor idx 122 (0-based) = 123 (1-based), which agrees.
CROP_FILE = "egfr_d3_crop_renum.pdb"
TGT_CROP_PATCH  = "120,121,122,123,124,125,126"   # canonical 430..436, H433 centred
TGT_CROP_ANCHOR = "123"                            # the anchor alone

ARMS = [
    ("tight_h433",    "DE", TGT_H433),   # the main play: winner's binder, narrowed target site
    ("tight_h433_DD", "DD", TGT_H433),   # chemistry sweep
    ("tight_h433_EE", "EE", TGT_H433),
    ("tight_h433_D1", "D1", TGT_H433),   # one carboxylate only -- is the second doing anything?
    ("dom4_pair",     "DE", TGT_DOM4),   # the only route past 5.41x
]

# (name, binder key, target indices) -- these use CROP_FILE instead of the full ECD
CROP_ARMS = [
    ("crop_patch",  "DE", TGT_CROP_PATCH),    # A/B vs tight_h433: crop is the only change
    ("crop_anchor", "DE", TGT_CROP_ANCHOR),   # declare only H433, to test whether MORE hint hurts
]

TEMPLATE = """\
# {name} -- generated by bin/make_tight_specs.py, see that file for the rationale.
# Binder band 56..90 (the only band that produced a switch). Pinned at {pins}.
# Target site (positional in egfr_ecd_6aru_renum.pdb): {tgt}
# Run with EXTRA_ARGS="--diffusion_batch_size 1" and DETACH=1.
entities:
  - protein:
      id: B
      sequence: {seq}
      binding_types:
        binding: {pins}
      secondary_structure:
        loop: {loops}

  - file:
      path: egfr_ecd_6aru_renum.pdb
      include:
        - chain:
            id: A
      binding_types:
        - chain:
            id: A
            binding: {tgt}
      structure_groups: "all"
"""


def spec(name, binder_key, tgt):
    seq, pins, loops = BINDER[binder_key]
    return TEMPLATE.format(name=name, seq=seq, pins=pins, loops=loops, tgt=tgt)


def selftest():
    for name, bk, tgt in ARMS + CROP_ARMS:
        seq, pins, loops = BINDER[bk]
        pinned = [int(x) for x in pins.split(",")]
        toks = [t.strip() for t in seq.split(",")]

        # letters uppercase, and the pin count matches the declared indices
        letters = [t for t in toks if t.isalpha()]
        assert all(t.isupper() for t in letters), f"{name}: lowercase pin becomes UNK"
        assert len(letters) == len(pinned), f"{name}: {len(letters)} letters vs {len(pinned)} indices"

        # no a..b range may precede the last pinned residue
        last_letter = max(i for i, t in enumerate(toks) if t.isalpha())
        for t in toks[:last_letter]:
            assert ".." not in t, f"{name}: range {t} precedes a pinned residue"

        # the pinned index is the running count of designed residues + 1, per letter
        n, idx = 0, []
        for t in toks:
            if t.isalpha():
                n += 1; idx.append(n)
            elif ".." in t:
                n += int(t.split("..")[0])
            else:
                n += int(t)
        assert idx == pinned, f"{name}: pins land at {idx}, spec declares {pinned}"

        # no loop window may cover a pinned residue (data.py:1953)
        for w in loops.split(","):
            a, b = (int(x) for x in w.split(".."))
            for p in pinned:
                assert not (a <= p <= b), f"{name}: loop {w} covers pinned residue {p}"

        # binder binding_types must be exactly the pins (data.py:1949)
        assert sorted(pinned) == sorted(int(x) for x in pins.split(",")), f"{name}: binding != pins"

        # H370 must be gone from every target site
        assert "343" not in tgt.split(","), f"{name}: H370 (343) is dead weight, drop it"

        # the declared length band must match BAND -- the header text quotes it, and a
        # comment that drifts from the spec is how the 123-160 band got rerun by mistake
        lo = sum(int(t.split("..")[0]) if ".." in t else (1 if t.isalpha() else int(t))
                 for t in toks)
        hi = sum(int(t.split("..")[1]) if ".." in t else (1 if t.isalpha() else int(t))
                 for t in toks)
        assert (lo, hi) == BAND, f"{name}: band is {(lo, hi)}, expected {BAND}"

        # every declared target index must be CENTRED on its anchor histidine, not merely
        # contain it. TGT_DOM4 was 554..558 (anchor 557 at position 4 of 5) until the
        # 2026-10-03 audit caught it; dom4_pair had already run off-centre.
        idx = sorted(int(x) for x in tgt.split(","))
        runs_ = []
        for i in idx:
            if runs_ and i == runs_[-1][-1] + 1: runs_[-1].append(i)
            else: runs_.append([i])
        for blk in runs_:
            if len(blk) >= 3:
                mid = blk[len(blk) // 2]
                assert mid in ANCHORS, (
                    f"{name}: target block {blk[0]}..{blk[-1]} is centred on {mid}, "
                    f"which is not an anchor histidine {sorted(ANCHORS)}")
    print(f"selftest OK: {len(ARMS)} arms")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest(); sys.exit(0)
    selftest()
    out = Path("targets/egfr")
    for name, bk, tgt in ARMS:
        p = out / f"boltzgen_egfr_{name}.yaml"
        p.write_text(spec(name, bk, tgt))
        print(f"wrote {p}")
    for name, bk, tgt in CROP_ARMS:
        p = out / f"boltzgen_egfr_{name}.yaml"
        p.write_text(spec(name, bk, tgt).replace("egfr_ecd_6aru_renum.pdb", CROP_FILE))
        print(f"wrote {p}  (cropped target)")
