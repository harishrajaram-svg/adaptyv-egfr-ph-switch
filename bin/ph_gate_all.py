#!/usr/bin/env python3
"""Run the RANKED OBJECTIVE on every design the project has generated.

The pH ratio is what the organizers rank first, and computing it costs 0.44 s of CPU
per design -- PROPKA twice (binder present, binder deleted in place from the SAME
coordinates) and a closed-form linkage. On 1,588 designs that is ~12 min on one core
and ~1 min on 18. Every narrowing this project has done -- generator iptm, then H433
contact distance, then an ipSAE screen -- was a cheap proxy standing in for a number
that was always affordable on the whole pool. This script removes the excuse.

It reports EVERY target histidine found, not a chosen one. Per site:
    ratio = K(6.5)/K(7.4),  K(pH) = (1+10^(pKa_bound-pH)) / (1+10^(pKa_free-pH))
plus the best single site and the product over all sites that help (ratio > 1).

Free leg: the binder chain is deleted IN PLACE, in reverse index order, from the same
file -- so every other coordinate is held fixed and the difference is the binder and
nothing else. An apo-reference free leg produced 10 false switches on 2026-10-03.
"""
import glob, os, sys, tempfile, subprocess, csv, json
from concurrent.futures import ProcessPoolExecutor
import gemmi

PH_LO, PH_HI = 6.5, 7.4
# UniProt name per residue number, per target family (chain-B HIS fingerprint).
CROP = {48:'H358',60:'H370',73:'H383',108:'H418',123:'H433'}
D3   = {24:'H358',36:'H370',49:'H383',84:'H418',99:'H433'}
ECD  = {23:'H47',121:'H145',159:'H183',209:'H233',280:'H304',334:'H358',346:'H370',
        359:'H383',394:'H418',409:'H433',483:'H507',535:'H559',560:'H584',566:'H590',
        591:'H615',594:'H618',597:'H621'}
# PROBLEM 2, registered 2026-10-08. Mature TNF-alpha 6-157, the gated canonical trimer; each
# of the three protomers carries the same three histidines. H73 is the one BinderBench names as
# a hotspot for this target, and it sits at the INTER-PROTOMER seam -- which is why the free leg
# must delete only the binder and keep the sibling protomers. See ph_gate_multisite.orient_multi.
TNF  = {15:'H15', 73:'H73', 78:'H78'}

def link(free, bound, ph_lo=None, ph_hi=None):
    """Thermodynamic linkage ratio K(ph_lo)/K(ph_hi).

    PARAMETERISED 2026-10-08. The pH pair used to be read from this module's globals while
    bin/ph_gate_multisite.py kept its OWN PH_LO/PH_HI copy. Setting the multisite copy for
    problem 2 (6.0/7.4) would have left this function silently computing problem 1's
    6.5/7.4 -- a wrong answer with no error. Callers now pass the pair; the defaults are
    problem 1's and are unchanged, so every existing result reproduces byte for byte.
    """
    lo = PH_LO if ph_lo is None else ph_lo
    hi = PH_HI if ph_hi is None else ph_hi
    K = lambda ph: (1 + 10**(bound-ph)) / (1 + 10**(free-ph))
    return K(lo)/K(hi)

def pkas(st, wd, tag, chain='B'):
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

def family(B):
    """Identify the target construct by its HISTIDINE FINGERPRINT, not its length.

    The length test used to be `len(B)==609`, which silently excluded the 621-residue
    ECD construct the scoring faa actually uses -- same numbering, 12 more residues --
    so every full-ECD refold was dropped without a word (32 designs, found 9:45 AM).
    The HIS tuple is the identity; match on that and the length is irrelevant.
    """
    hs=tuple(r.seqid.num for r in B if r.name=='HIS')
    if hs==tuple(sorted(CROP)): return CROP
    if hs==tuple(sorted(D3)):   return D3
    if hs==tuple(sorted(ECD)):  return ECD
    if hs==tuple(sorted(TNF)):  return TNF
    return None

def orient(model):
    """Return (target_chain, family) by TESTING BOTH CHAIN ORDERS.

    Chain order is NOT consistent across this project. 29 of 69 ESMFold2 run directories
    put the TARGET in chain A and the binder in chain B -- every `reval_*` and
    `rimA01_r02_s*` run among them, which are revalidation runs for SUBMITTED designs --
    while the BoltzGen generator poses and all of today's runs put the binder first.
    Assuming A=binder skipped 150 poses, and on a construct that happened to pass the
    fingerprint on the wrong chain it would have measured the binder's own histidines
    silently. The histidine fingerprint decides which chain is the target.
    """
    C={c.name:c for c in model}
    for tgt in ('B','A'):
        if tgt in C and len(C)>1:
            f=family(C[tgt])
            if f is not None: return tgt,f
    return None,None

def one(cif):
    try:
        st=gemmi.read_structure(cif); st.setup_entities()
        st.remove_ligands_and_waters(); st.setup_entities()
        tgt,names=orient(st[0])
        if names is None: return None
        with tempfile.TemporaryDirectory() as wd:
            bound=pkas(st,wd,'cpx',tgt)
            m=st[0]                              # same file, BINDER removed in place
            for i in range(len(m)-1,-1,-1):
                if m[i].name!=tgt: del m[i]
            st.setup_entities()
            free=pkas(st,wd,'apo',tgt)
        sites={}
        for num,nm in names.items():
            f,b=free.get(num),bound.get(num)
            if f is None or b is None: continue
            sites[nm]={'free':round(f,2),'bound':round(b,2),'ratio':round(link(f,b),3)}
        if not sites: return None
        # `product` MUST compose over ALL sites, not only the ones that help.
        # Until 2026-10-04 it multiplied only the ratio>1.0 sites, which is selection on the
        # outcome: a site whose ratio is <1 is a site where acid genuinely weakens binding and
        # it belongs in the product. 1,547 of 1,584 designs had an overstated value, by up to
        # 4.77x, and the error manufactured the project's only claim of beating the 7.943x
        # single-site maximum -- rimA01/d3_rimA_50 was reported at 9.10x (H370 6.056 x H383
        # 1.503) while its own H433 reads 0.723. Its true all-site product is 6.581x, and with
        # the fix NO design in the pool exceeds the single-site ceiling.
        # Sites PROPKA never moved return exactly 1.000 and are harmless in either form.
        helpful=[v['ratio'] for v in sites.values() if v['ratio']>1.0]
        prod=1.0
        for v in sites.values(): prod*=v['ratio']
        best=max(sites.items(), key=lambda kv: kv[1]['ratio'])
        parts=cif.split('/')
        return {'arm':parts[2],'file':os.path.basename(cif)[:-4],'path':cif,
                'target_chain':tgt,
                'best_site':best[0],'best_ratio':best[1]['ratio'],
                'product':round(prod,3),           # all sites, see above
                'product_helpful_only':round(__import__('math').prod(helpful) if helpful else 1.0,3),
                'n_helpful':len(helpful),'sites':sites}
    except Exception as e:
        return None

if __name__=='__main__':
    files=sorted(glob.glob('runs/*/*/intermediate_designs_inverse_folded/refold_cif/*.cif'))
    files=[f for f in files if '1alu-smoketest' not in f]
    print(f'gating {len(files)} designs on every target histidine', flush=True)
    rows=[]
    with ProcessPoolExecutor(max_workers=16) as ex:
        for i,r in enumerate(ex.map(one, files, chunksize=8),1):
            if r: rows.append(r)
            if i%200==0: print(f'  {i}/{len(files)}', flush=True)
    json.dump(rows, open('analysis/01-egfr/ph_gate_all.json','w'))
    print(f'\n{len(rows)} designs gated, {len(files)-len(rows)} skipped (unknown target family)')
    rows.sort(key=lambda r:-r['best_ratio'])
    print(f'\n=== TOP 40 BY BEST SINGLE-SITE RATIO ===')
    print(f"{'ratio':>7} {'site':>6} {'prod':>7} {'arm':<20} design")
    for r in rows[:40]:
        print(f"{r['best_ratio']:7.2f} {r['best_site']:>6} {r['product']:7.2f} {r['arm']:<20} {r['file']}")
    import collections
    print(f'\n=== which histidine is the best site, across the pool ===')
    for s,n in collections.Counter(r['best_site'] for r in rows).most_common():
        print(f'  {s:>6}  {n}')
    print(f'\ndesigns with best_ratio >= 4.0: {sum(1 for r in rows if r["best_ratio"]>=4.0)}')
    print(f'designs with best_ratio >= 5.0: {sum(1 for r in rows if r["best_ratio"]>=5.0)}')
    print(f'designs with >=2 helpful sites: {sum(1 for r in rows if r["n_helpful"]>=2)}')
