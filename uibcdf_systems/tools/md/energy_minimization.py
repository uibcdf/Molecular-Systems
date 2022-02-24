
def energy_minimization(item, platform_name='CUDA', verbose=True):

    from openmm import LangevinIntegrator, Platform, Context, LocalEnergyMinimizer_minimize
    from openmm import unit

    # Integrator.

    integrator = LangevinIntegrator(0*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)

    # Platform.

    platform = Platform.getPlatformByName(platform_name)

    # Context.

    context = Context(item.system, integrator, platform)
    context.setPositions(item.coordinates)

    # Minimization.

    if verbose==True:
        energy = context.getState(getEnergy=True).getPotentialEnergy()
        print('Potential energy before minimization: {}'.format(energy))

    LocalEnergyMinimizer_minimize(context)

    if verbose==True:
        energy = context.getState(getEnergy=True).getPotentialEnergy()
        print('Potential energy after minimization: {}'.format(energy))

    item.coordinates = context.getState(getPositions=True).getPositions(asNumpy=True)

    pass
