"""TNF-alpha epitope census, computed not recalled.
Contact = heavy atom within 4.5 A of a heavy atom of a non-TNF polymer chain.
Chain ID and numbering offset are both solved by matching observed residues to P01375
BY SEQID (gap-safe), requiring >=0.90 identity over >=40 observed residues."""
import gemmi, sys, json
from collections import defaultdict
HUMAN="VRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTHTISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL"
CUT=4.5
def one(r):
    i=gemmi.find_tabulated_residue(r.name)
    return i.one_letter_code.upper() if i and i.is_amino_acid() else None
def best_offset(ch):
    """k such that uniprot = seqid + k. Returns (k, identity, n_compared)."""
    obs=[(r.seqid.num,one(r)) for r in ch if one(r)]
    if len(obs)<40: return (None,0,0)
    best=(None,0,0)
    for k in range(40,120):
        m=t=0
        for num,c in obs:
            u=num+k-77
            if 0<=u<len(HUMAN): t+=1; m+=(c==HUMAN[u])
        if t>=40 and m/t>best[1]: best=(k,m/t,t)
    return best
out={}
for pdb in sys.argv[1:]:
    st=gemmi.read_structure(pdb+".cif"); st.setup_entities(); st.remove_ligands_and_waters()
    model=st[0]
    tnf={}; other=[]
    for ch in model:
        k,acc,n=best_offset(ch)
        if acc>=0.90: tnf[ch.name]=k
        elif sum(1 for r in ch if one(r))>=15: other.append(ch.name)
    if not tnf or not other:
        print(f"{pdb}: SKIP (tnf {len(tnf)}, partner {len(other)})"); continue
    ns=gemmi.NeighborSearch(model,st,5.0).populate() if False else gemmi.NeighborSearch(st,5.0).populate()
    # contacts per (partner chain -> tnf chain -> set of uniprot resnums)
    per=defaultdict(lambda: defaultdict(set)); allhits={}
    for ch in model:
        if ch.name not in tnf: continue
        k=tnf[ch.name]
        for res in ch:
            aa=one(res)
            if aa is None: continue
            u=res.seqid.num+k
            for atom in res:
                if atom.element==gemmi.Element('H'): continue
                for m in ns.find_atoms(atom.pos,'\0',radius=CUT):
                    cra=m.to_cra(model)
                    if cra.chain.name in other and cra.atom.element!=gemmi.Element('H') \
                       and cra.atom.pos.dist(atom.pos)<=CUT:
                        per[cra.chain.name][ch.name].add(u); allhits[u]=aa
    uni=sorted(allhits)
    # per-partner two-protomer split (take the partner chain with the most contacts)
    rich=max(per, key=lambda p: sum(len(v) for v in per[p].values()))
    split={c:sorted(per[rich][c]) for c in per[rich]}
    out[pdb]={"offsets":tnf,"partners":other,"epitope":[f"{allhits[u]}{u}" for u in uni],
              "n":len(uni),"split_example":{"partner":rich,"per_tnf_chain":split}}
    print(f"\n### {pdb}  TNF chains {sorted(tnf)} (uniprot=seqid+{sorted(set(tnf.values()))}), partners {other}")
    print(f"  union epitope, {len(uni)} residues:\n    "+" ".join(f"{allhits[u]}{u}" for u in uni))
    print(f"  one partner chain ({rich}) splits across protomers: "+
          "; ".join(f"{c}:{len(v)}" for c,v in sorted(split.items(), key=lambda x:-len(x[1]))))
json.dump(out,open("epitopes.json","w"),indent=1)
