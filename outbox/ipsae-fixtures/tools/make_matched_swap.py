#!/usr/bin/env python3
"""Build case 12: a TRUE matched chain-order swap of case 01.

On 2026-10-05 the reviewer pointed out that case 06 is an altogether different complex which
merely happens to list its chains the other way round -- it is not one complex and its PAE
re-ordered against itself. Correct. Case 06 shows the pipeline handles a binder-first file,
but it has no target-first counterpart, so it cannot show that ipSAE_min is
INVARIANT to chain order -- which is the property the submission depends on,
since 29 of 69 run directories place the target in chain A and 40 do not.

This builds the matched counterpart from case 01 by permuting the file, not by
picking a second structure:
  * residues are re-emitted B-first then A-first;
  * chain labels are exchanged (old B -> A, old A -> B), in both label_asym_id
    and auth_asym_id;
  * the PAE matrix and the pLDDT vector are permuted with the SAME index
    permutation, so every residue keeps its own PAE row/column.
Identical coordinates, identical PAE content, opposite chain order. ipSAE_min
must come out bit-identical; anything else is an order dependence.
"""
import json, shutil, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
SRC = HERE / "cases" / "01_barnase_barstar_positive_regression"
DST = HERE / "cases" / "12_matched_chain_order_swap"
STEM_IN = "pos_barnase_barstar_sample_0"
STEM_OUT = "swap_barnase_barstar_sample_0"

def parse_atom_site(text):
    hdr, rows, pre, post, inloop, done = [], [], [], [], False, False
    for ln in text.splitlines():
        if ln.startswith("_atom_site."):
            hdr.append(ln.strip()); inloop = True; continue
        if inloop and not done:
            if ln.startswith("#") or (ln.strip() and ln.startswith("_")):
                done = True; post.append(ln); continue
            if ln.strip():
                rows.append(ln.split()); continue
            done = True; post.append(ln); continue
        (post if done else pre).append(ln)
    return pre, hdr, rows, post

def main():
    cif = (SRC / f"{STEM_IN}.cif").read_text()
    pre, hdr, rows, post = parse_atom_site(cif)
    LA, LS = hdr.index("_atom_site.label_asym_id"), hdr.index("_atom_site.label_seq_id")
    AA = hdr.index("_atom_site.auth_asym_id")
    AS = hdr.index("_atom_site.auth_seq_id")
    AID = hdr.index("_atom_site.id")

    # residue order as the file presents it -- this is the PAE index order
    order, idx_of = [], {}
    for r in rows:
        k = (r[LA], r[LS])
        if k not in idx_of:
            idx_of[k] = len(order); order.append(k)
    chains = []
    for c, _ in order:
        if not chains or chains[-1] != c: chains.append(c)
    if len(chains) != 2:
        sys.exit(f"expected 2 chains in file order, got {chains}")
    first, second = chains
    SWAP = {first: second, second: first}

    # permutation: emit `second` chain's residues first, then `first`
    new_order = [k for k in order if k[0] == second] + [k for k in order if k[0] == first]
    perm = [idx_of[k] for k in new_order]            # perm[new] = old

    # rewrite atom rows in the new residue order, exchanging chain labels
    by_res = {}
    for r in rows: by_res.setdefault((r[LA], r[LS]), []).append(r)
    out_rows, serial = [], 0
    for k in new_order:
        for r in by_res[k]:
            r = list(r); serial += 1
            r[LA] = SWAP[r[LA]]; r[AA] = SWAP[r[AA]]; r[AID] = str(serial)
            out_rows.append(r)

    DST.mkdir(parents=True, exist_ok=True)
    widths = [max(len(r[i]) for r in out_rows) for i in range(len(hdr))]
    body = "\n".join(" ".join(f.ljust(w) for f, w in zip(r, widths)).rstrip() for r in out_rows)
    (DST / f"{STEM_OUT}.cif").write_text(
        "\n".join(pre) + "\n" + "\n".join(hdr) + "\n" + body + "\n" + "\n".join(post) + "\n")

    # permute PAE + pLDDT with the same permutation
    d = json.loads((SRC / f"{STEM_IN}_ipsae.json").read_text())
    pae = d["pae"]
    n = len(pae)
    if n != len(order):
        sys.exit(f"PAE is {n}x{n} but the file has {len(order)} residues -- refusing")
    d["pae"] = [[pae[perm[i]][perm[j]] for j in range(n)] for i in range(n)]
    if isinstance(d.get("plddt"), list) and len(d["plddt"]) == n:
        d["plddt"] = [d["plddt"][perm[i]] for i in range(n)]
    (DST / f"{STEM_OUT}_ipsae.json").write_text(json.dumps(d))

    print(f"wrote {DST}")
    print(f"  file order {first}({sum(1 for k in order if k[0]==first)}) -> "
          f"{second}({sum(1 for k in order if k[0]==second)}) became "
          f"{SWAP[second]} -> {SWAP[first]}")
    print(f"  PAE permuted {n}x{n}; residue 0 of the new file was residue {perm[0]} of the old")

if __name__ == "__main__":
    main()
