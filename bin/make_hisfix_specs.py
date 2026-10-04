#!/usr/bin/env python3
"""Generate the Mechanism A pinned-residue BoltzGen specs for problem 1's final arms.

Eight specs from one template, because hand-writing eight near-identical YAMLs is how the
CIF converter ended up with four divergent copies and one silent bug.

DESIGN (2026-10-03). The measured failure was that 106 of 140 Mechanism A designs (76%)
had NO binder histidine at the interface -- pH lived in the filter, never in the design
objective. These specs move it. Factorial:

    labels OFF vs ON   : arms 1-2 vs 3-4   (identical sequence spec; labels are the only change)
    one His vs two      : arms 3-4 vs 5-6
    His vs Gln          : arms 3-4 vs 7-8   (Gln is near-isosteric, NO titratable proton)
    domII vs domIII     : within every pair

WHAT THE LABELS CAN AND CANNOT SAY (BoltzGen 247b9bb, verified in source):
  * binding_types on a DESIGNED protein entity is the FLAT dict form, and it is UNARY --
    "this residue is interfacial", partner unspecified. The model's only pairwise channel
    (contact_conditioning, featurizer.py:788-804) is hard-coded to UNSPECIFIED with no YAML
    route, so "the His must touch E455" is NOT expressible. Writing the file:-style
    list-of-chain form here parses and does NOTHING (schema.py:870-900).
  * secondary_structure indexes the SAMPLED chain, so the pinned residue's index must be
    deterministic -- every a..b range has to come AFTER it. That is why the sequence specs
    are "60, Z, 59..99" and not "60..80, Z, 59..79": with a leading range no index can
    address the pinned residue at all.
  * Target binding indices are POSITIONAL within the chain as loaded, not author resseq.
    egfr_ecd_6aru_renum.pdb is contiguous 1..609, so positional == renum. This is only safe
    BECAUSE the file was renumbered; check contiguity before reusing on another target.
  * Case is load-bearing and silent: lowercase 'h' becomes UNK, not HIS (const.py:163-191).
"""
import sys
from pathlib import Path

SITES = {
    "domII":  dict(binding="292,293,294", canon="E319/E320/D321", dom="II",
                   note="three CONTIGUOUS acidic residues, the most coupled protons on the "
                        "receptor; 20.5 A sequon clearance, 3/3 conserved, burial 23. "
                        "Caveat: 16.9 A from the Y246-D563 tether, and domain II rearranges "
                        "between the tethered and extended states."),
    "domIII": dict(binding="428,431,433", canon="E455/D458/D460", dom="III",
                   note="20.6 A sequon clearance, 3/3 conserved, and 41.9 A from the tether "
                        "-- the farthest of any clear site, so insensitive to the "
                        "tethered-vs-extended question. Domain III is also the drug-relevant "
                        "surface (cetuximab epitope, ligand binding). The 2026-10-01 scan "
                        "wrongly concluded domain III had no usable acidic surface."),
    "neutral": dict(binding="260,267,270", canon="F287/K294/R297", dom="II",
                   note="THE NEGATIVE CONTROL for the carboxylate-proximity requirement. "
                        "That requirement was called necessary off 6/6 switches having a "
                        "carboxylate within 4.0 A and 0/96 non-switches switching -- but "
                        "every Mechanism A design we ever made was aimed AT an acidic site, "
                        "so the claim had no negative arm. Caveat stated plainly: a truly "
                        "carboxylate-free patch does NOT exist on this surface. The EGFR ECD "
                        "is acid-dense -- median nearest-carboxylate distance for a "
                        "non-acidic residue is 5.6 A and the maximum anywhere is 17.2 A, and "
                        "every position beyond 9 A is either a structural cysteine or a "
                        "glycosylation sequon. This patch is the best available: worst-case "
                        "7.5 A to any carboxylate, 1.9x the 4.0 A the mechanism requires, "
                        "conserved, burial 19-45, sequon clearance >20 A."),
}
# (label, sequence spec, pinned indices in the SAMPLED chain, loop windows)
# Loop windows must EXCLUDE the pinned indices. BoltzGen validates
# (src/boltzgen/data/data.py:1953, DesignInfo.is_valid) that every residue carrying a
# secondary-structure label is also DESIGNED -- and a pinned letter is fixed, so it is not.
# A window spanning the pinned residue dies with "Misspecified design info. There were
# residues that have a secondary structure type specified but are not set to be designed."
# Flanking the pinned residue with loop on both sides achieves the same thing: it sits in a
# flexible, solvent-exposed stretch without itself carrying a label.
PINS = {
    "1His": ("60, H, 59..99",        "61",    "56..60,62..66"),
    "2His": ("45, H, 28, H, 45..85", "46,75", "41..45,47..51,70..74,76..80"),
    "1Gln": ("60, Q, 59..99",        "61",    "56..60,62..66"),
}
ARMS = [("domII", "1His", False), ("domIII", "1His", False),      # labels OFF (control)
        ("domII", "1His", True),  ("domIII", "1His", True),
        ("domII", "2His", True),  ("domIII", "2His", True),
        ("domII", "1Gln", True),  ("domIII", "1Gln", True),       # non-titratable control
        ("neutral", "1His", True), ("neutral", "2His", True)]     # no-carboxylate control


def seq_len(spec):
    """Total binder length range implied by a sequence spec. Fixed letters count 1 each."""
    lo = hi = 0
    for part in (p.strip() for p in spec.split(",")):
        if ".." in part:
            a, b = part.split(".."); lo += int(a); hi += int(b)
        elif part.isdigit():
            lo += int(part); hi += int(part)
        else:
            lo += len(part); hi += len(part)
    return lo, hi


def render(site, pin, labels):
    s, p = SITES[site], PINS[pin]
    spec, idx, loops = p
    lo, hi = seq_len(spec)
    assert (lo, hi) == (120, 160), f"{pin}: length {lo}..{hi}, must match the 120..160 baseline"
    res = "HISTIDINE" if "His" in pin else "GLUTAMINE (non-titratable control)"
    extra = ""
    if labels:
        extra = (f"\n      binding_types:\n        binding: {idx}"
                 f"\n      secondary_structure:\n        loop: {loops}")
    return f"""# MECHANISM A, {res} PINNED AT GENERATION -- {site} site, {pin}, labels {'ON' if labels else 'OFF'}.
#
# Across 140 Mechanism A designs on two sites, 106 (76%) contained NO binder histidine at
# the interface and exactly ONE was a real switch. The limiting factor was never the site or
# the mechanism: nothing told the generator to supply a titratable group. This spec does.
#
# Site: {s['canon']} canonical (renum {s['binding']}), domain {s['dom']}.
#   {s['note']}
#
# Pinned residue at sampled index {idx}. Length {lo}..{hi}, identical to the baseline, so the
# pinned residue (and the labels, where present) are the only changes.
# Labels {'ON: the pinned residue is declared interfacial and placed in a loop. The interface' if labels else 'OFF: this is the attribution control -- same sequence spec as its labels-ON'}
# {'label is UNARY (partner unspecified); a pairwise "contacts E455" term does not exist.' if labels else 'twin, differing only by the absence of binding_types and secondary_structure.'}
#
# DIGITS request that many designed residues, LETTERS are fixed identities (schema.py
# 743-760). No wildcard exists: 'X' pins UNK, lowercase silently becomes UNK.
# Run with EXTRA_ARGS="--diffusion_batch_size 1".
# Numbering: canonical = renum + 27, mature = renum + 3, verified residue-by-residue.
entities:
  - protein:
      id: B
      sequence: {spec}{extra}

  - file:
      path: egfr_ecd_6aru_renum.pdb
      include:
        - chain:
            id: A
      binding_types:
        - chain:
            id: A
            binding: {s['binding']}
      structure_groups: "all"
"""


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        for k, (spec, idx, loops) in PINS.items():
            assert seq_len(spec) == (120, 160), (k, seq_len(spec))
            # Mirrors make_g532mimic_specs.py:128-131. Absent these three, a spec with a
            # lowercase pin (silently UNK) or a range before a pinned residue (index
            # re-sampled, so no label can address it) passed the self-test unchallenged.
            assert all(c.isupper() for c in spec if c.isalpha()), f"{k}: lowercase becomes UNK"
            pinned = [int(x) for x in idx.split(",")]
            for part in spec.split(",")[: 2 * len(pinned)]:
                assert ".." not in part, f"{k}: a range precedes a pinned residue -- index not addressable"
            # a loop window must not cover a pinned letter (data.py:1953)
            for w in loops.split(","):
                a, b = (int(x) for x in w.split(".."))
                for q in pinned:
                    assert not (a <= q <= b), f"{k}: loop {w} covers pinned residue {q}"
        assert seq_len("60, H, 59..99") == (120, 160)
        assert seq_len("45, H, 28, H, 45..85") == (120, 160)
        assert len(ARMS) == 10 and len({a[:2] + (a[2],) for a in ARMS}) == 10
        # labels-OFF arms must be byte-identical to their labels-ON twin minus the two blocks
        on, off = render("domII", "1His", True), render("domII", "1His", False)
        # match the YAML BLOCKS, not the substrings -- the labels-OFF comment text
        # legitimately mentions binding_types and secondary_structure by name.
        assert "\n      binding_types:\n        binding: 61" in on
        assert "\n      secondary_structure:\n        loop:" in on
        assert "\n      binding_types:\n        binding: 61" not in off
        assert "\n      secondary_structure:" not in off
        assert off.split("sequence:")[1].split("\n")[0] == on.split("sequence:")[1].split("\n")[0]
        def expand(w):
            out=set()
            for part in w.split(","):
                a,b=part.split(".."); out |= set(range(int(a), int(b)+1))
            return out
        for k,(spec, idx, loops) in PINS.items():
            pinned={int(i) for i in idx.split(",")}
            overlap = pinned & expand(loops)
            assert not overlap, f"{k}: loop window covers pinned index {overlap} -- BoltzGen rejects SS on a non-designed residue"
        print("self-test OK: all specs 120..160, 10 distinct arms, control differs only by "
              "labels, no loop window covers a pinned residue")
        sys.exit()
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "targets/egfr")
    for site, pin, labels in ARMS:
        name = f"boltzgen_egfr_hf_{site}_{pin}{'' if labels else '_ctrl'}.yaml"
        (out / name).write_text(render(site, pin, labels))
        print(f"  {name}")
