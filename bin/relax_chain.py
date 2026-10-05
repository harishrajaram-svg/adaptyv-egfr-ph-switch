#!/usr/bin/env python3
"""Side-chain relaxation of a single chain extracted from a bound complex.

This is the "relaxed" leg of the reviewer's sensitivity request:
  "the limitation is the retained bound conformation: SIDE-CHAIN RELAXATION, water
   penetration, conformational populations and protonation coupling in the actual free
   protein are not represented."

It sits between the two legs already computed, and that is the point of it:

  deletion    partner's atoms removed, everything else frozen in the bound pose
  RELAXED     partner removed, then side chains minimised with the BACKBONE RESTRAINED
  apo fold    the chain folded alone by ESMFold2, a different backbone entirely

Restraining the backbone is deliberate. Released, this would become a slow refold and
would no longer isolate the one effect being tested. So the relaxed leg answers exactly
one question: how much of the deletion estimate is an artefact of side chains that never
got to move?

Implicit solvent (GBn2) rather than vacuum: in vacuum, surface polar side chains collapse
onto the protein to satisfy their own electrostatics, which is precisely the kind of
burial change a pKa calculation is sensitive to, and would manufacture the effect.

    relax_chain.py <complex.cif> <keep_chain> <out.pdb> [--steps 500]
"""
import sys


def relax(cif, keep_chain, out_pdb, steps=500, restraint_k=10.0):
    import gemmi
    from io import StringIO
    import openmm, openmm.app as app
    from openmm import unit
    from pdbfixer import PDBFixer

    st = gemmi.read_structure(str(cif))
    st.setup_entities(); st.remove_ligands_and_waters(); st.remove_alternative_conformations()
    keep = [ch for ch in st[0] if ch.name == keep_chain]
    if not keep:
        raise SystemExit(f"{cif}: no chain {keep_chain} (have {[c.name for c in st[0]]})")
    out = gemmi.Structure(); out.add_model(gemmi.Model('1'))
    c = gemmi.Chain('A')
    for r in keep[0]:
        c.add_residue(r)
    out[0].add_chain(c); out.setup_entities()
    tmp = StringIO(out.make_pdb_string())

    fixer = PDBFixer(pdbfile=tmp)
    fixer.findMissingResidues()
    fixer.missingResidues = {}          # do NOT build loops; relax what is there
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    ff = app.ForceField('amber14/protein.ff14SB.xml', 'implicit/gbn2.xml')
    system = ff.createSystem(fixer.topology, nonbondedMethod=app.NoCutoff,
                             constraints=app.HBonds)
    # harmonic restraint on backbone heavy atoms
    force = openmm.CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    force.addGlobalParameter('k', restraint_k * unit.kilocalories_per_mole / unit.angstroms**2)
    for p in ('x0', 'y0', 'z0'):
        force.addPerParticleParameter(p)
    n = 0
    for atom in fixer.topology.atoms():
        if atom.name in ('N', 'CA', 'C', 'O'):
            pos = fixer.positions[atom.index]
            force.addParticle(atom.index, [pos.x, pos.y, pos.z])
            n += 1
    system.addForce(force)

    integrator = openmm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                                 0.002 * unit.picoseconds)
    sim = app.Simulation(fixer.topology, system, integrator,
                         openmm.Platform.getPlatformByName('CPU'))
    sim.context.setPositions(fixer.positions)
    sim.minimizeEnergy(maxIterations=steps)
    pos = sim.context.getState(getPositions=True).getPositions()
    with open(out_pdb, 'w') as fh:
        app.PDBFile.writeFile(fixer.topology, pos, fh, keepIds=True)
    return n


if __name__ == '__main__':
    a = [x for x in sys.argv[1:] if not x.startswith('--')]
    steps = int(sys.argv[sys.argv.index('--steps') + 1]) if '--steps' in sys.argv else 500
    n = relax(a[0], a[1], a[2], steps=steps)
    print(f"relaxed {a[0]} chain {a[1]} -> {a[2]} ({n} backbone atoms restrained)")
