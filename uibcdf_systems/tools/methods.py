def get_potential_energy(item, coordinates=None, platform_name='CUDA'):

    from openmm import LangevinIntegrator, Platform, Context
    from openmm import unit
    import numpy as np

    integrator = LangevinIntegrator(0.0*unit.kelvin, 0.0/unit.picoseconds, 2.0*unit.femtoseconds)

    platform = Platform.getPlatformByName(platform_name)

    context = Context(item.system, integrator, platform)

    if coordinates is None:
        context.setPositions(item.coordinates)
    else:
        context.setPositions(coordinates)

    if item.box is not None:
        context.setPeriodicBoxVectors(item.box[0], item.box[1], item.box[2])

    state = context.getState(getEnergy=True)
    potential_energy = state.getPotentialEnergy()

    return potential_energy

def get_probability_density(item, range=None, n_bins=None, temperature=None):

    from openmm import unit

    if type(item)==unit.Quantity:

        from numpy import histogram

        length_unit = item.unit
        region = [range[0]._value, range[1]._value]
        probability_density, bin_limites = histogram(item._value, bins=n_bins, range=region, density=True)
        coordinates = (bin_limites[1:]+bin_limites[:-1])*0.50
        delta_x = (bin_limites[1:]-bin_limites[:-1]).mean()

        probability_density = probability_density/length_unit
        bin_limites = bin_limites*length_unit
        coordinates = coordinates*length_unit
        delta_x = delta_x*length_unit

    elif hasattr(item, 'system'):

        from numpy import histogram_bin_edges, exp, zeros

        Z = 0.0
        kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
        kBT = kB * temperature

        length_unit = range[0].unit
        lim_inf = range[0]._value
        lim_sup = range[1]._value

        bin_limites = histogram_bin_edges([lim_inf, lim_sup], n_bins)*length_unit
        coordinates = (bin_limites[1:]+bin_limites[:-1])*0.50
        delta_x = (bin_limites[1:]-bin_limites[:-1]).mean()

        coordinates_aux = zeros((coordinates.shape[0],3))*length_unit
        coordinates_aux[:,0] = coordinates[:]
        boltzmann_weights = exp(-item.potential(coordinates_aux)/kBT)
        Z = boltzmann_weights.sum()*delta_x
        probability_density = boltzmann_weights/Z

    return probability_density, coordinates, bin_limites, delta_x


def get_partition_function(item, range=None, n_bins=None, temperature=None):

    if hasattr(item, 'system'):

        from numpy import histogram_bin_edges, exp

        Z = 0.0
        kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
        kBT = kB * temperature

        length_unit = range[0].unit
        lim_inf = range[0]._value
        lim_sup = range[1]._value

        bin_limites = histogram_bin_edges([lim_inf, lim_sup], n_bins)*length_unit
        coordinates = (bin_limites[1:]+bin_limites[:-1])*0.50*length_unit
        delta_x = (bin_limites[1:]-bin_limites[:-1]).mean()*length_unit

        Z = (exp(-item.potential(coordinates)/kBT)).sum()*delta_x

    return Z

