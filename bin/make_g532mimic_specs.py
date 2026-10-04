#!/usr/bin/env python3
"""Arms that copy G532's actual mechanism, plus rim-placement arms. 2026-10-03.

WHY. By tonight three independent legs say the same thing: good binding means burial, and
burial destroys the pH switch (mechanism B on H433; pinned-histidine mechanism A at noise
level with controls matching; and 0/59 switches in the best-bound acid-site designs, which
carry the LARGEST pKa suppressions we have measured, -0.9 to -2.6).

But G532 gets 8-13x anyway (Liu et al. 2022, PMID 36458200). Three features of it we never
combined:
  1. the CARBOXYLATES are on the binder and it titrates the TARGET's histidines
     -- we pinned a histidine on the binder instead;
  2. it engages BOTH target histidines, H433 and H370, which sit 7.7 A apart (side chain,
     measured on egfr_ecd_6aru_renum) and are rigid across seven structures -- we only ever
     targeted H433. This matters because our binder-side 2-His arms came out WORSE than
     1-His (p=0.027): coupled protons on a floppy designed binder interfere, whereas two
     protons on the rigid target are pre-organised;
  3. it reaches in with a CDR loop -- a shallow rim contact, not a buried core.

The rim arms test (3) directly by labelling the pinned residue `not_binding` while its
sequence neighbours are `binding`. If burial is what kills the switch, that is the single
change that should rescue it.

Numbering: renum = canonical - 27 = mature - 3, verified residue-by-residue.
  canonical H370 -> renum 343, canonical H433 -> renum 406 (both HIS, confirmed)
  domII acidic cluster E319/E320/D321 -> renum 292,293,294 (CLEAR, 20.5 A sequon clearance)

Constraints that are load-bearing and were learned the hard way today:
  * every a..b range must come AFTER the pinned residue or its index is re-sampled and no
    label can address it;
  * a loop window must NOT cover a pinned residue -- BoltzGen validates that every
    SS-labelled residue is also designed, and a pinned letter is not (data.py:1953);
  * letters must be UPPERCASE; lowercase silently becomes UNK.
"""
import sys
from pathlib import Path

TGT_HIS  = "343,406"      # canonical H370 + H433 -- for binder-side CARBOXYLATE arms
TGT_ACID = "292,293,294"  # canonical E319/E320/D321 -- for binder-side HISTIDINE arms

# name -> (sequence, pinned indices, loop windows, target, rim?)
ARMS = {
 "g532mimic":       ("45, D, 28, E, 45..85", "46,75", "41..45,47..51,70..74,76..80", TGT_HIS,  False),
 "g532mimic_short": ("20, D, 14, E, 20..54", "21,36", "16..20,22..26,31..35,37..41", TGT_HIS,  False),
 "g532mimic_tiny":  ("12, D, 8, E, 8..28",   "13,22", "9..12,14..17,19..21,23..26",  TGT_HIS,  False),
 # The two `rim_*` arms are REMOVED. They labelled designed binder residues as `binding`,
 # which BoltzGen rejects: "Only target residues can have a binding type specified since
 # this feature indicates where the design should bind" (data.py:1949). Rim placement is
 # NOT expressible through binding_types. The three g532mimic arms test the same
 # hypothesis by SIZE instead -- 120-160 / 56-90 / 30-50 residues is a burial
 # dose-response, since a 30-mer cannot bury a residue the way a 160-mer can.
 "tinyHis":         ("12, H, 17..37",        "13",    "9..12,14..17",                TGT_ACID, False),
}
EXPECT_LEN = {"g532mimic": (120, 160), "g532mimic_short": (56, 90), "g532mimic_tiny": (30, 50),
              "tinyHis": (30, 50)}


def seq_len(spec):
    lo = hi = 0
    for part in (p.strip() for p in spec.split(",")):
        if ".." in part:
            a, b = part.split(".."); lo += int(a); hi += int(b)
        elif part.isdigit():
            lo += int(part); hi += int(part)
        else:
            lo += len(part); hi += len(part)
    return lo, hi


def expand(w):
    out = set()
    for part in w.split(","):
        a, b = part.split(".."); out |= set(range(int(a), int(b) + 1))
    return out


def render(name):
    spec, idx, loops, tgt, rim = ARMS[name]
    pinned = [int(i) for i in idx.split(",")]
    lo, hi = seq_len(spec)
    mech = "B (binder carboxylate titrates a TARGET histidine)" if tgt == TGT_HIS \
           else "A (binder histidine paired with a TARGET carboxylate)"
    if rim:
        # neighbours interfacial, the pinned residue explicitly NOT -- rim, not core
        lo_n, hi_n = pinned[0] - 6, pinned[0] + 6
        bind = f"{lo_n}..{pinned[0]-1},{pinned[0]+1}..{hi_n}"
        bt = f"\n      binding_types:\n        binding: {bind}\n        not_binding: {idx}"
        rimnote = ("\n# RIM ARM: the pinned residue is labelled not_binding while its neighbours are\n"
                   "# binding, which places the titratable group at the interface RIM rather than in\n"
                   "# its core. This is the direct test of the burial mechanism.")
    else:
        bt = f"\n      binding_types:\n        binding: {idx}"
        rimnote = ""
    return f"""# {name} -- mechanism {mech}
# Length {lo}..{hi}. Pinned at sampled index/indices {idx}, which is deterministic because
# every a..b range sits after them.{rimnote}
# Target hotspots (renum): {tgt}  [{'canonical H370+H433' if tgt == TGT_HIS else 'canonical E319/E320/D321'}]
# Run with EXTRA_ARGS="--diffusion_batch_size 1".
entities:
  - protein:
      id: B
      sequence: {spec}{bt}
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


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        for n, (spec, idx, loops, tgt, rim) in ARMS.items():
            assert seq_len(spec) == EXPECT_LEN[n], (n, seq_len(spec), EXPECT_LEN[n])
            pinned = {int(i) for i in idx.split(",")}
            bad = pinned & expand(loops)
            assert not bad, f"{n}: loop window covers pinned {bad} -- BoltzGen rejects SS on a fixed residue"
            assert spec.split(",")[0].strip().isdigit(), f"{n}: first block must be a fixed count"
            # everything up to and including the LAST pinned letter must be range-free;
            # the trailing block is allowed (and expected) to be a range.
            for part in spec.split(",")[: 2 * len(pinned)]:
                assert ".." not in part, f"{n}: a range precedes a pinned residue -- index not addressable"
            assert all(c.isupper() for c in spec if c.isalpha()), f"{n}: lowercase becomes UNK"
        # binding_types may ONLY name the pinned (non-designed) residue, never a
        # designed one -- BoltzGen rejects the latter outright.
        for n, (spec, idx, loops, tgt, rim) in ARMS.items():
            r = render(n)
            assert "not_binding" not in r, f"{n}: not_binding is unsupported here"
            bl = r.split("binding_types:")[1].split("binding:")[1].split("\n")[0].strip()
            assert bl == idx, f"{n}: binder binding_types must be exactly the pinned indices, got {bl!r}"
        print(f"self-test OK: {len(ARMS)} arms, lengths exact, no loop covers a pinned residue, "
              "no range precedes one, all uppercase")
        sys.exit()
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "targets/egfr")
    for n in ARMS:
        (out / f"boltzgen_egfr_{n}.yaml").write_text(render(n))
        print(f"  boltzgen_egfr_{n}.yaml")
