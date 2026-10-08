#!/usr/bin/env python3
"""pae_interface_min / _mean -- the inter-chain PAE metrics, as a ranker.

WHY THIS EXISTS, and why it is a sibling of ipsae_min.py rather than an edit to it.
analysis/02-tnf/calibrate_external.py measures every ranking metric this project could
use against 150 TNF-alpha designs with wet-lab labels. On OUR OWN validated arm,
ESMFold2-Full, at the 1-binder-to-trimer construct we actually submit:

    pae_interface_min   AUC 0.901        <-- best on our arm
    ipsae_max           AUC 0.899
    ipsae_min           AUC 0.822        <-- what we rank on today
    iptm_pae            AUC 0.833
    plddt_binder        AUC 0.579        (a filter, not a ranker)

So the switch is pre-registered on external measured data, not on our own designs:
rank on pae_interface_min, keep ipSAE_min as the reported sibling. ipsae_min.py is the
most load-bearing module in the repo and is pinned by 12 fixtures plus a fail-closed
suite; adding a metric inside it would put that suite at risk for no gain.

LOWER IS BETTER. These are PAE values in angstroms, so they are negated before any
ranking that assumes higher-is-better. calibrate_external.py does that explicitly.

FAIL CLOSED. A score is returned only when the chain labelling read off the OUTPUT
structure accounts for every row of the PAE matrix. The project's standing trap list
says chain labels change through the folding wrapper, so they are never assumed from
the input.

    pae_interface.py <pae_json> <structure.cif>
    pae_interface.py --dir <dir>      # every *_ipsae.json next to its *.cif
    pae_interface.py --selftest
"""
import glob
import json
import os
import sys
from pathlib import Path


def chain_ids(struct_path):
    """Per-residue chain label, in PAE row order: model order, polymer residues only."""
    import gemmi
    st = gemmi.read_structure(str(struct_path))
    st.setup_entities()
    out = []
    for ch in st[0]:
        poly = ch.get_polymer()
        n = len(poly) if len(poly) else len(ch)
        out.extend([ch.name] * n)
    return out


def interface_pae(pae, chains):
    """(min, mean, n_pairs) over residue pairs whose two chains differ.

    Both directions are included: PAE is not symmetric, and the protocol's ipSAE
    definition takes the min over directions for the same reason.
    """
    n = len(pae)
    if len(chains) != n:
        return None
    vals = [pae[i][j] for i in range(n) for j in range(n) if chains[i] != chains[j]]
    if not vals:
        return None
    return min(vals), sum(vals) / len(vals), len(vals)


def score(pae_file, struct):
    try:
        d = json.load(open(pae_file))
    except Exception as e:
        print(f"  REFUSE {Path(pae_file).name}: unreadable PAE ({type(e).__name__})")
        return None
    pae = d.get("pae")
    if not pae or not isinstance(pae, list) or not isinstance(pae[0], list):
        print(f"  REFUSE {Path(pae_file).name}: no square 'pae' matrix")
        return None
    try:
        chains = chain_ids(struct)
    except Exception as e:
        print(f"  REFUSE {Path(struct).name}: unreadable structure ({type(e).__name__})")
        return None
    r = interface_pae(pae, chains)
    if r is None:
        print(f"  REFUSE {Path(struct).name}: {len(chains)} chain labels for "
              f"{len(pae)} PAE rows, or a single chain -- no interface to score")
        return None
    return r


def main():
    args = sys.argv[1:]
    if args and args[0] == "--dir":
        pairs = []
        for j in sorted(glob.glob(os.path.join(args[1], "**", "*_ipsae.json"), recursive=True)):
            cif = Path(j.replace("_ipsae.json", ".cif"))
            if cif.exists():
                pairs.append((Path(j), cif))
        if not pairs:
            sys.exit(f"no *_ipsae.json + .cif pairs under {args[1]}")
    elif len(args) >= 2:
        pairs = [(Path(args[0]), Path(args[1]))]
    else:
        sys.exit(__doc__)

    print(f"{'structure':44}{'pae_if_min':>11}{'pae_if_mean':>12}{'pairs':>9}")
    ok = 0
    for pae, cif in pairs:
        r = score(pae, cif)
        if r is None:
            continue
        mn, mu, npair = r
        ok += 1
        print(f"{cif.name[:43]:44}{mn:>11.3f}{mu:>12.3f}{npair:>9}")
    print(f"\nscored {ok} of {len(pairs)}.  LOWER IS BETTER.")
    return 0 if ok else 1


def selftest():
    # A 4-residue, 2-chain complex. Inter-chain block holds 1.0 and 9.0; the
    # intra-chain diagonal blocks hold 0.1, which must never be reached.
    pae = [
        [0.1, 0.1, 1.0, 9.0],
        [0.1, 0.1, 5.0, 5.0],
        [1.0, 5.0, 0.1, 0.1],
        [9.0, 5.0, 0.1, 0.1],
    ]
    chains = ["A", "A", "B", "B"]
    mn, mu, n = interface_pae(pae, chains)
    assert n == 8, n                      # 2x2 block, both directions
    assert mn == 1.0, mn
    assert abs(mu - (1.0 + 9.0 + 5.0 + 5.0) * 2 / 8) < 1e-12, mu
    # MUTATION 1: counting intra-chain pairs too would reach the 0.1 diagonal.
    allpairs = [pae[i][j] for i in range(4) for j in range(4) if i != j]
    assert min(allpairs) == 0.1, "fixture is wrong"
    assert mn != min(allpairs), "intra-chain pairs must be excluded, and are not"
    # MUTATION 2: a one-chain structure has no interface and must REFUSE, not return 0.
    assert interface_pae(pae, ["A", "A", "A", "A"]) is None
    # MUTATION 3: a chain list that does not account for every PAE row must REFUSE.
    assert interface_pae(pae, ["A", "A", "B"]) is None
    assert interface_pae(pae, ["A", "A", "B", "B", "B"]) is None
    # asymmetry is preserved: swapping one direction changes the min
    asym = [row[:] for row in pae]
    asym[0][2] = 0.4
    assert interface_pae(asym, chains)[0] == 0.4
    # three chains: every unlike pair counts, which is the 1:trimer construct
    tri = [[0.1] * 3 for _ in range(3)]
    tri[0][1] = tri[1][0] = 2.0
    tri[0][2] = tri[2][0] = 3.0
    tri[1][2] = tri[2][1] = 4.0
    assert interface_pae(tri, ["A", "B", "C"]) == (2.0, 3.0, 6)
    print("pae_interface.py --selftest PASS")
    return 0


if __name__ == "__main__":
    sys.exit(selftest() if "--selftest" in sys.argv else main())
