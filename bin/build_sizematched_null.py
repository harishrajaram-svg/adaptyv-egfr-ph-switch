#!/usr/bin/env python3
"""Size-matched shuffled negative for the pae_interface control band (item 2f).

The 4 existing negatives (targets/validation-tnf/neg_shuffled_tnfr2_*) shuffle the 164-residue
TNFR2 ectodomain, so they score at 302,382 cross-chain residue pairs. The 35 problem-2 designs
are 76 or 84 residues and score at 219,486 / 227,022 pairs. pae_interface_mean is an average over
those pairs, so a 164-mer negative is not a yardstick for a 76-mer design -- the control band
committed in p2_control_band_2026-10-08.tsv cannot be compared to the designs as it stands.

This builds one negative at the designs' own length AND own composition: chains A/B/C are the
same 157-mer TNF-alpha protomers read from the gated fixture (not retyped), chain D is a
composition-matched shuffle of a named design, same convention and seed as FIXTURES.json
records for the TNFR2 shuffles.

The design is named on the command line and must be chosen BEFORE scoring. The one used is
p2trimer02_L76_s0: the first L76 row of the committed candidate TSV, picked by file order so the
choice cannot be read off the metric being calibrated.
"""
import csv
import random
import sys

FIXTURE = "targets/validation-tnf/pos_tnf_tnfr2.faa"
CANDIDATES = "outbox/02-tnf-candidates.tsv"
SEED = 20261006          # FIXTURES.json negatives.seed


def read_faa(path):
    chains, name = [], None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            chains.append([name, ""])
        elif line:
            chains[-1][1] += line
    return chains


def trimer(path=FIXTURE):
    """The 3 TNF-alpha protomers from the gated fixture. Refuses anything else."""
    chains = read_faa(path)
    abc = chains[:3]
    if [c[0] for c in abc] != ["protein|A", "protein|B", "protein|C"]:
        sys.exit(f"REFUSE: {path} chains are not A/B/C: {[c[0] for c in chains]}")
    seqs = [c[1] for c in abc]
    if len(set(seqs)) != 1 or len(seqs[0]) != 157:
        sys.exit(f"REFUSE: A/B/C are not 3 identical 157-mers: {[len(s) for s in seqs]}")
    return seqs[0]


def design_seq(name, path=CANDIDATES):
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["design"] == name:
                seq, n = row["sequence"], int(row["length"])
                if len(seq) != n:
                    sys.exit(f"REFUSE: {name} length column {n} != sequence {len(seq)}")
                return seq
    sys.exit(f"REFUSE: {name} not in {path}")


def shuffled(seq, seed=SEED):
    out = list(seq)
    random.Random(seed).shuffle(out)
    out = "".join(out)
    if out == seq:
        sys.exit("REFUSE: shuffle returned the input unchanged")
    if sorted(out) != sorted(seq):
        sys.exit("REFUSE: shuffle did not preserve composition")
    return out


def build(name, out_path):
    tnf, null = trimer(), shuffled(design_seq(name))
    with open(out_path, "w") as fh:
        for cid in "ABC":
            fh.write(f">protein|{cid}\n{tnf}\n")
        fh.write(f">protein|D\n{null}\n")
    total = 3 * len(tnf) + len(null)
    pairs = total ** 2 - (3 * len(tnf) ** 2 + len(null) ** 2)
    print(f"wrote {out_path}: 3x{len(tnf)} + {len(null)} = {total} residues, "
          f"{pairs} cross-chain pairs (shuffle of {name}, seed {SEED})")


def selftest():
    tnf = trimer()
    seq = design_seq("p2trimer02_L76_s0")
    sh = shuffled(seq)
    assert len(sh) == len(seq) == 76, (len(sh), len(seq))
    assert sorted(sh) == sorted(seq), "composition not preserved"
    assert sh != seq, "order not destroyed"
    assert shuffled(seq) == sh, "not reproducible at a fixed seed"
    # the whole point: cross-chain pair count must equal what the L76 designs scored at
    total = 3 * len(tnf) + len(sh)
    pairs = total ** 2 - (3 * len(tnf) ** 2 + len(sh) ** 2)
    assert pairs == 219486, f"pair count {pairs} != the 219486 the L76 designs scored at"
    # MUTATION TEST: an 84-mer null must NOT land on the L76 pair count
    n84 = 3 * len(tnf) + 84
    assert n84 ** 2 - (3 * len(tnf) ** 2 + 84 ** 2) == 227022
    print(f"  ok  trimer read from the fixture: 3x{len(tnf)}")
    print(f"  ok  chain D is a composition-matched, reproducible shuffle of a 76-mer")
    print(f"  ok  cross-chain pairs {pairs} == the L76 designs' own pair count")
    print(f"  ok  mutation test: an 84-mer would be 227022, not {pairs}")
    print("\nself-tests passed: 4")


if __name__ == "__main__":
    if len(sys.argv) == 1 or "--selftest" in sys.argv:
        selftest()
    else:
        build(sys.argv[1], sys.argv[2])
