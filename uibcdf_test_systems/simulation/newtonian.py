
def newtonian(system, friction=None,
              initial_positions=None, initial_velocities=None, integration_timestep=None,
              saving_timestep=None, total_time=None, platform_name='CPU', verbose=True):


    from simtk.openmm import LangevinIntegrator, Platform, Context
    from simtk import unit
    import numpy as np

    # System parameters.
    n_particles = system.getNumParticles()

    # Integration parameters.

    steps_per_cicle = int(round(saving_timestep/integration_timestep))
    n_steps = int(round(total_time/integration_timestep))
    n_cicles = int(round(n_steps/steps_per_cicle))

    # Integrator.

    temperature = 0.0*unit.kelvin
    if friction is None:
        friction = 0.0/unit.picoseconds
    integrator = LangevinIntegrator(temperature, friction, integration_timestep)

    # Platform.

    platform = Platform.getPlatformByName(platform_name)

    # Context.

    context = Context(system, integrator, platform)
    context.setPositions(initial_positions)
    context.setVelocities(initial_velocities)

    # Reporter arrays: time, position, velocity, kinetic_energy, potential_energy

    time = np.zeros([n_cicles], np.float32) * unit.picoseconds
    position = np.zeros([n_cicles, n_particles, 3], np.float32) * unit.nanometers
    velocity = np.zeros([n_cicles, n_particles, 3], np.float32) * unit.nanometers/unit.picosecond
    kinetic_energy = np.zeros([n_cicles, n_particles, 3], np.float32) * unit.kilocalories_per_mole
    potential_energy = np.zeros([n_cicles, n_particles, 3], np.float32) * unit.kilocalories_per_mole

    # Initial context in reporters

    state = context.getState(getPositions=True, getVelocities=True, getEnergy=True)
    time[0] = state.getTime()
    position[0] = state.getPositions()
    velocity[0] = state.getVelocities()
    kinetic_energy[0] = state.getKineticEnergy()
    potential_energy[0] = state.getPotentialEnergy()

    # Integration loop saving every cicle steps

    for ii in range(1, n_cicles):
        context.getIntegrator().step(steps_per_cicle)
        state = context.getState(getPositions=True, getVelocities=True, getEnergy=True)
        time[ii] = state.getTime()
        position[ii] = state.getPositions()
        velocity[ii] = state.getVelocities()
        kinetic_energy[ii] = state.getKineticEnergy()
        potential_energy[ii] = state.getPotentialEnergy()

    return time, position, velocity, kinetic_energy, potential_energy

