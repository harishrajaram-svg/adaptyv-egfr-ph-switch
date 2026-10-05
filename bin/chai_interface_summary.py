#!/usr/bin/env python3
"""Does Chai-1 form an interface where ESMFold2 did not? The independent check.

Reviewer, 2026-10-05: "If feasible, run a bounded independent structure-prediction check
on the failed positives and a diverse finalist subset; ESMFold2-Fast alone is not an
independent validation."

WHAT CAN AND CANNOT BE COMPARED. Chai-1 emits no residue-level PAE -- its score file
carries aggregate_score, ptm, iptm, per_chain_pair_iptm and clash flags only -- so
ipSAE_min cannot be computed on its structures and there is no like-for-like SCORE
comparison to be had. The comparison is therefore GEOMETRIC: does an interface form, how
large is it, and (for rAC1, which has a solved complex) is it the right one. Chai's own
iptm is reported alongside but this project has already measured ipTM as near-chance on
this target, so it is metadata, not evidence.

Chai-1 also uses ESM embeddings by default, so "independent" here means an independent
architecture, weights and training procedure -- not independence from ESM sequence
featurisation. Disabling them would degrade the model without achieving that, so the
validated configuration was used.

    chai_interface_summary.py [--cut 5.0]
"""
import glob, json, os, sys
from collections import defaultdict
import gemmi


# Minimal .npz reader, so this script has no numpy dependency.
#
# A .npz is a zip of .npy members; each .npy is a short ASCII header dict followed by raw
# little-endian data. The values needed here are all shape-(1,) float32 or bool, so a
# struct.unpack of the first element is enough. Keeping numpy out means the published tree
# runs with gemmi alone, and the fixture bundle already showed what an avoidable import
# costs a reviewer.
def npz_scalars(path, keys):
    import struct, zipfile
    out = {}
    with zipfile.ZipFile(path) as z:
        for name in z.namelist():
            key = name[:-4] if name.endswith('.npy') else name
            if key not in keys:
                continue
            raw = z.read(name)
            if raw[:6] != b'\x93NUMPY':
                continue
            hlen = struct.unpack('<H', raw[8:10])[0] if raw[6] == 1 else \
                   struct.unpack('<I', raw[8:12])[0]
            hstart = 10 if raw[6] == 1 else 12
            hdr = raw[hstart:hstart + hlen].decode('latin1')
            body = raw[hstart + hlen:]
            if "'<f4'" in hdr or "'|f4'" in hdr:
                out[key] = struct.unpack('<f', body[:4])[0]
            elif "'|b1'" in hdr or "'<b1'" in hdr:
                out[key] = bool(body[0])
    return out

ROOT = 'runs/chai1/w4_indep'
CUT = 5.0


def chains(path):
    st = gemmi.read_structure(str(path))
    st.setup_entities(); st.remove_ligands_and_waters()
    return [(ch.name, [r for r in ch]) for ch in st[0] if len(ch) >= 10]


def n_contacts(a, b, cut):
    cell = cut
    grid = defaultdict(list)
    for rb in b:
        for at in rb:
            if at.element == gemmi.Element('H'):
                continue
            grid[(int(at.pos.x//cell), int(at.pos.y//cell), int(at.pos.z//cell))].append(at.pos)
    pairs = set()
    for i, ra in enumerate(a):
        for at in ra:
            if at.element == gemmi.Element('H'):
                continue
            kx, ky, kz = int(at.pos.x//cell), int(at.pos.y//cell), int(at.pos.z//cell)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for p in grid.get((kx+dx, ky+dy, kz+dz), ()):
                            if at.pos.dist(p) <= cut:
                                pairs.add(i); break
    return len(pairs)


def main():
    cut = CUT
    if '--cut' in sys.argv:
        cut = float(sys.argv[sys.argv.index('--cut') + 1])
    man_p = ('/private/tmp/claude-501/-Users-harish-code-context-directory/'
             '8a78545e-a035-42db-b995-8251eebdb653/scratchpad/chai_in/manifest.json')
    notes = {m['tag']: m['note'] for m in json.load(open(man_p))} if os.path.exists(man_p) else {}
    rows = []
    print(f"{'complex':<40}{'models':>7}{'iptm':>7}{'clash':>6}{'iface_res':>10}  note")
    for d in sorted(glob.glob(os.path.join(ROOT, '*'))):
        if not os.path.isdir(d):
            continue
        tag = os.path.basename(d)
        cifs = sorted(glob.glob(os.path.join(d, 'pred.model_idx_*.cif')))
        if not cifs:
            continue
        iptms, clashes, ifaces = [], 0, []
        for c in cifs:
            npz = c.replace('pred.', 'scores.').replace('.cif', '.npz')
            if os.path.exists(npz):
                z = npz_scalars(npz, {'iptm', 'has_inter_chain_clashes'})
                if 'iptm' in z:
                    iptms.append(z['iptm'])
                clashes += int(z.get('has_inter_chain_clashes', False))
            ch = chains(c)
            if len(ch) >= 2:
                ch.sort(key=lambda x: -len(x[1]))
                # interface = target residues contacting ANY other chain
                other = [r for nm, rs in ch[1:] for r in rs]
                ifaces.append(n_contacts(ch[0][1], other, cut))
        import statistics as _st
        med_iptm = float(_st.median(iptms)) if iptms else float('nan')
        med_if = int(_st.median(ifaces)) if ifaces else 0
        rows.append(dict(tag=tag, n_models=len(cifs), iptm_median=round(med_iptm, 4),
                         clash_models=clashes, iface_residues_median=med_if,
                         iface_residues_all=ifaces, note=notes.get(tag, '')))
        print(f"{tag[:39]:<40}{len(cifs):>7}{med_iptm:>7.3f}{clashes:>6}{med_if:>10}  "
              f"{notes.get(tag,'')[:46]}")
    json.dump(dict(cut=cut, complexes=rows),
              open('analysis/01-egfr/chai_interface_summary.json', 'w'), indent=1)
    print(f"\nwrote analysis/01-egfr/chai_interface_summary.json ({len(rows)} complexes)")
    print("NOTE: no ipSAE is computed on Chai structures -- Chai emits no residue-level PAE.")


if __name__ == '__main__':
    main()
