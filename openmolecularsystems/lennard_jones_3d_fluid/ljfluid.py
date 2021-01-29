from simtk import openmm
from simtk import unit
from simtk.openmm import app
import numpy as np


def subrandom_particle_positions(nparticles, box_vectors, ndim):
    """Generate a deterministic list of subrandom particle positions."""
    # Create positions array.
    positions = unit.Quantity(np.zeros([nparticles, 3], np.float32), unit.nanometers)

    # Generate Sobol' sequence.
    import sobol
    ivec = sobol.i4_sobol_generate(ndim, nparticles, 1)
    x = np.array(ivec, np.float32)
    for dim in range(ndim):
        l = box_vectors[dim][dim]
        positions[:, dim] = unit.Quantity((x[dim, :] - 0.5) * l / l.unit, l.unit)

    return positions


def distributeAtoms(boxsize=[34, 34, 34],
                     nparticles=1000,
                     reduced_density=0.05,
                     mass=39.9 * unit.amu,  # argon
                     sigma=3.4 * unit.angstrom,  # argon,
                     epsilon=0.238 * unit.kilocalories_per_mole,  # argon,
                     cutoff=None,
                     switch_width=3.4 * unit.angstrom,  # argon
                     dispersion_correction=True,
                     lattice=False,
                     charge=None,
                     **kwargs):
        # Determine Lennard-Jones cutoff.
        if cutoff is None:
            cutoff = 3.0 * sigma

        charge = 0.0 * unit.elementary_charge
        cutoff_type = openmm.NonbondedForce.CutoffPeriodic

        # Create an empty system object.
        system = openmm.System()

        # Periodic box vectors.
        a = unit.Quantity((boxsize[0] * unit.angstrom, 0 * unit.angstrom, 0 * unit.angstrom))
        b = unit.Quantity((0 * unit.angstrom, boxsize[1] * unit.angstrom, 0 * unit.angstrom))
        c = unit.Quantity((0 * unit.angstrom, 0 * unit.angstrom, boxsize[2] * unit.angstrom))
        system.setDefaultPeriodicBoxVectors(a, b, c)

        # Set up periodic nonbonded interactions with a cutoff.
        nb = openmm.NonbondedForce()
        nb.setNonbondedMethod(cutoff_type)
        nb.setCutoffDistance(cutoff)
        nb.setUseDispersionCorrection(dispersion_correction)

        nb.setUseSwitchingFunction(False)
        if (switch_width is not None):
            nb.setUseSwitchingFunction(True)
            nb.setSwitchingDistance(cutoff - switch_width)

        for particle_index in range(nparticles):
            system.addParticle(mass)
            nb.addParticle(charge, sigma, epsilon)

        positions = subrandom_particle_positions(nparticles, system.getDefaultPeriodicBoxVectors(), 2)
        # Add the nonbonded force.
        system.addForce(nb)

        # Add a restrining potential to keep atoms in z=0
        energy_expression = 'k * (z^2)'
        force = openmm.CustomExternalForce(energy_expression)
        force.addGlobalParameter('k', 100)
        for particle_index in range(nparticles):
            force.addParticle(particle_index, [])
        system.addForce(force)

        # Create topology.
        topology = app.Topology()
        element = app.Element.getBySymbol('Ar')
        chain = topology.addChain()
        for particle in range(system.getNumParticles()):
            residue = topology.addResidue('Ar', chain)
            topology.addAtom('Ar', element, residue)
        topology.setUnitCellDimensions(unit.Quantity(boxsize, unit.angstrom)) 
            
        # Simulate it
        from simtk.openmm import LangevinIntegrator, VerletIntegrator
        from simtk.openmm.app import Simulation, PDBReporter, StateDataReporter, PDBFile
        from simtk.unit import kelvin, picoseconds, picosecond, angstrom
        from sys import stdout
        from mdtraj.reporters import DCDReporter
        #from dcdreporter import DCDReporter
        nsteps = 10000
        freq = 1
        #integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
        integrator = VerletIntegrator(0.002 * picoseconds)
        simulation = Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)
        simulation.minimizeEnergy()
        simulation.reporters.append(DCDReporter('output.dcd', 1))
        simulation.reporters.append(StateDataReporter(stdout, 1000, potentialEnergy=True, totalEnergy=True, step=True, separator='   '))
        simulation.step(nsteps)

        state = simulation.context.getState(getPositions=True)
        finalpos = state.getPositions(asNumpy=True).value_in_unit(angstrom)

        with open('topology.pdb', 'w') as f:
            PDBFile.writeFile(topology, positions, f)

        from htmd.molecule.molecule import Molecule
        mol = Molecule('topology.pdb')
        mol.read('output.dcd')

        return finalpos, mol, system, simulation


if __name__ == '__main__':
    finalpos = distributeAtoms(reduced_density=1)
