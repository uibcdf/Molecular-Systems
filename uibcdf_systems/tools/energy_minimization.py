import openmm as mm
import openmm.unit as unit

def energy_minimization(item, platform_name='CUDA', verbose=False):

    # Integrator.

    integrator = mm.LangevinIntegrator(0*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)

    # Platform.

    platform = mm.Platform.getPlatformByName(platform_name)

    # Context.

    context = mm.Context(item.system, integrator, platform)
    context.setPositions(item.coordinates)

    # Minimization.

    if verbose==True:
        energy = context.getState(getEnergy=True).getPotentialEnergy()
        print('Potential energy before minimization: {}'.format(energy))

    LocalEnergyMinimizer_minimize(context)

    if verbose==True:
        energy = context.getState(getEnergy=True).getPotentialEnergy()
        print('Potential energy after minimization: {}'.format(energy))

    coordinates = context.getState(getPositions=True).getPositions(asNumpy=True)
    item.set_coordinates(coordinates)

    pass

