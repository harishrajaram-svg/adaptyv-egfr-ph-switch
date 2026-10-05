import os, sys
import gemmi, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fetch
from collections import defaultdict
HUMAN="VRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTHTISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL"
CUT=4.5
def one(r):
    i=gemmi.find_tabulated_residue(r.name)
    return i.one_letter_code.upper() if i and i.is_amino_acid() else None
def best_offset(ch):
    obs=[(r.seqid.num,one(r)) for r in ch if one(r)]
    if len(obs)<40: return (None,0)
    best=(None,0)
    for k in range(-20,130):
        m=t=0
        for num,c in obs:
            u=num+k-77
            if 0<=u<len(HUMAN): t+=1; m+=(c==HUMAN[u])
        if t>=40 and m/t>best[1]: best=(k,m/t)
    return best
LABEL = {}   # labels.json was a scratchpad convenience and is not an input
res_out={}
for pdb in sys.argv[1:]:
    st=gemmi.read_structure(fetch.cif(pdb)); st.setup_entities(); st.remove_ligands_and_waters()
    model=st[0]; tnf={}; other={}
    for ch in model:
        k,acc=best_offset(ch)
        if acc>=0.90: tnf[ch.name]=k
        elif sum(1 for r in ch if one(r))>=15: other[ch.name]=sum(1 for r in ch if one(r))
    if not tnf or not other: print(f"{pdb}: SKIP tnf={len(tnf)} other={len(other)}"); continue
    ns=gemmi.NeighborSearch(st,5.0).populate()
    per=defaultdict(lambda: defaultdict(set)); aa={}
    for ch in model:
        if ch.name not in tnf: continue
        k=tnf[ch.name]
        for r in ch:
            c=one(r)
            if c is None: continue
            u=r.seqid.num+k
            for at in r:
                if at.element==gemmi.Element('H'): continue
                for m in ns.find_atoms(at.pos,'\0',radius=CUT):
                    cra=m.to_cra(model)
                    if cra.chain.name in other and cra.atom.element!=gemmi.Element('H') and cra.atom.pos.dist(at.pos)<=CUT:
                        per[cra.chain.name][ch.name].add(u); aa[u]=c
    # group partner chains that belong to one binding unit: take each partner chain separately,
    # report the single richest one and its two-protomer split
    print(f"\n### {pdb}  offset +{sorted(set(tnf.values()))[0]}  TNF:{len(tnf)} chains  partners:{sorted(other)}")
    ranked=sorted(per, key=lambda p:-sum(len(v) for v in per[p].values()))
    for p in ranked[:4]:
        tot=set().union(*per[p].values())
        split="; ".join(f"{c}({len(v)})" for c,v in sorted(per[p].items(), key=lambda x:-len(x[1])))
        print(f"  partner {p} [{other[p]}aa] -> {len(tot)} TNF residues across protomers {split}")
        print(f"      {' '.join(f'{aa[u]}{u}' for u in sorted(tot))}")
    res_out[pdb]={p:{"total":sorted(set().union(*per[p].values())),
                     "per_protomer":{c:sorted(v) for c,v in per[p].items()}} for p in ranked}
# write next to the script, not into whatever directory it was launched from (s26:
# a relative output path means the artifact lands somewhere different in a clean clone)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "per_partner.json")
json.dump(res_out, open(out, "w"), indent=1)
print(f"\nwrote {out}")
