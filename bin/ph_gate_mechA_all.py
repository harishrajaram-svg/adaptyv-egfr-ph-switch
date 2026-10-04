#!/usr/bin/env python3
"""Mechanism A across every design the project has generated.

Mechanism B puts the titratable group on the TARGET (EGFR's own histidines) and the
stabilising carboxylate on the binder. Mechanism A inverts it: the histidine is on the
BINDER, where we control its environment, and the negative charge is a target Asp/Glu.
Same linkage, opposite locus of control:

    pKa_bound > pKa_free  ->  protonation favours binding  ->  binds harder in acid

Conventions are taken from bin/ph_gate_mechA.py unchanged (PH 6.5/7.4, RATIO_BAR 1.20,
binder = chain A, target = chain B) so numbers from the two are comparable.

The free leg deletes the TARGET in place, in reverse index order, from the same file --
the mirror of what ph_gate_all.py does to the binder, and for the same reason: holding
every other coordinate fixed makes the difference the partner chain and nothing else.
Collecting wanted residues and re-adding them instead dangles references into the
chain's vector while it empties (the knockout_control.py bug, 2026-10-04).

Also reports the nearest TARGET carboxylate to each binder histidine ring. That is
mechanism A's geometry, the analogue of the H433 contact distance for mechanism B, and
the distance bins measured for B (canonical H-bond 2.6-3.2 A switching 62%) are the
prior worth testing here -- not assuming it carries over.
"""
import glob, os, sys, tempfile, subprocess, json
from concurrent.futures import ProcessPoolExecutor
import gemmi

PH_LO, PH_HI = 6.5, 7.4
RATIO_BAR = 1.20
TIP = {"ASP": ("OD1","OD2"), "GLU": ("OE1","OE2")}

def link(free, bound):
    K = lambda ph: (1 + 10**(bound-ph)) / (1 + 10**(free-ph))
    return K(PH_LO)/K(PH_HI)

def pkas(st, wd, tag, chain):
    st.write_pdb(os.path.join(wd, tag+'.pdb'))
    subprocess.run([sys.executable,'-m','propka',tag+'.pdb'],capture_output=True,cwd=wd)
    f=os.path.join(wd,tag+'.pka'); out={}
    if os.path.exists(f):
        for ln in open(f):
            q=ln.split()
            if len(q)>3 and q[0]=='HIS' and q[2]==chain:
                try: out[int(q[1])]=float(q[3])
                except ValueError: pass
    return out

def nearest_target_acid(A_res, B):
    """Closest target carboxylate oxygen to this binder histidine's ring nitrogens."""
    ring=[a.pos for a in A_res if a.name in ("ND1","NE2")]
    if not ring: return None, None
    best, who = 999.0, None
    for r in B:
        if r.name not in TIP: continue
        for a in r:
            if a.name in TIP[r.name]:
                for rn in ring:
                    d=a.pos.dist(rn)
                    if d<best: best, who = d, f"{r.name}{r.seqid.num}"
    return best, who

def one(cif):
    try:
        st=gemmi.read_structure(cif); st.setup_entities()
        st.remove_ligands_and_waters(); st.setup_entities()
        C={c.name:c for c in st[0]}
        A,B=C.get('A'),C.get('B')
        if A is None or B is None: return None
        his=[r for r in A if r.name=='HIS']
        if not his: return None                      # no binder histidine -> mechanism A impossible
        geom={r.seqid.num: nearest_target_acid(r,B) for r in his}
        with tempfile.TemporaryDirectory() as wd:
            bound=pkas(st,wd,'cpx','A')
            m=st[0]                                  # same file, TARGET removed in place
            for i in range(len(m)-1,-1,-1):
                if m[i].name!='A': del m[i]
            st.setup_entities()
            free=pkas(st,wd,'apo','A')
        sites={}
        for r in his:
            n=r.seqid.num
            f,b=free.get(n),bound.get(n)
            if f is None or b is None: continue
            d,who=geom[n]
            sites[f'HIS{n}']={'free':round(f,2),'bound':round(b,2),
                              'ratio':round(link(f,b),3),
                              'acid_dist':None if d is None else round(d,2),'acid':who}
        if not sites: return None
        best=max(sites.items(), key=lambda kv: kv[1]['ratio'])
        return {'arm':cif.split('/')[2],'file':os.path.basename(cif)[:-4],'path':cif,
                'n_his':len(sites),'best_his':best[0],'best_ratio':best[1]['ratio'],
                'best_acid':best[1]['acid'],'best_dist':best[1]['acid_dist'],'sites':sites}
    except Exception:
        return None

if __name__=='__main__':
    files=[f for f in sorted(glob.glob('runs/*/*/intermediate_designs_inverse_folded/refold_cif/*.cif'))
           if '1alu-smoketest' not in f]
    print(f'mechanism A over {len(files)} designs', flush=True)
    rows=[]
    with ProcessPoolExecutor(max_workers=16) as ex:
        for i,r in enumerate(ex.map(one, files, chunksize=8),1):
            if r: rows.append(r)
            if i%200==0: print(f'  {i}/{len(files)}', flush=True)
    json.dump(rows, open('analysis/01-egfr/ph_gate_mechA_all.json','w'))
    print(f'\n{len(rows)}/{len(files)} designs carry >=1 binder histidine and gated')
    rows.sort(key=lambda r:-r['best_ratio'])
    print(f'\n=== TOP 40 BY MECHANISM-A RATIO ===')
    print(f"{'ratio':>7} {'his':>8} {'acid':>9} {'dist':>6}  {'arm':<20} design")
    for r in rows[:40]:
        d='   -  ' if r['best_dist'] is None else f"{r['best_dist']:6.2f}"
        print(f"{r['best_ratio']:7.2f} {r['best_his']:>8} {str(r['best_acid']):>9} {d}  {r['arm']:<20} {r['file']}")
    for bar in (1.20, 2.0, 3.0, 4.0, 5.0):
        print(f"\ndesigns with mechanism-A ratio >= {bar}: {sum(1 for r in rows if r['best_ratio']>=bar)}")
    # does the mechanism-B distance prior carry over?
    print('\n=== switch rate by binder-HIS / target-acid distance (the B prior was 62% at 2.6-3.2 A) ===')
    BINS=[(0,2.4),(2.4,2.6),(2.6,3.2),(3.2,4.0),(4.0,6.0),(6.0,1e9)]
    for lo,hi in BINS:
        sub=[r for r in rows if r['best_dist'] is not None and lo<=r['best_dist']<hi]
        if not sub: continue
        sw=sum(1 for r in sub if r['best_ratio']>=RATIO_BAR)
        print(f'  {lo:4.1f}-{hi if hi<100 else 999:5.1f} A   n={len(sub):4d}   ratio>=1.20: {sw:4d} = {100*sw/len(sub):5.1f}%')
