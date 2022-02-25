from ._molecular_system import MolecularSystem
import openmm.unit as unit
import openmm as mm
import openmm.app as app
import numpy as np
from sympy import symbols

class HarmonicWell(MolecularSystem):

    """Particles in an harmonic well potential

    Test system with particles in an harmonic well potential.

    .. math::

        \\frac{1}{2} k\\left(x^2 + y^2 + z^2\\right)

    Attributes
    ----------
    n_particles
        Number of particles
    mass
        Mass of particles
    system
        Openmm system
    potential_expression
        External potential expression as a sympy function.
    potential_parameters
        Dictionary with the values of the parameters of the potential.

    Methods
    -------
    potential
        Potential evaluation at certain coordinates.

    """

    def __init__(self, n_particles=1, mass=32*unit.amu,
                 k=10.0*unit.kilocalorie/(unit.mole*unit.nanometers**2),
                 coordinates=None):

        """Creating a new instance of HarmonicWell

        A new test system is returned with the openmm system of particles in an external harmonic
        well potential:

        Parameters
        ----------

        n_particles: int
            Number of particles in the system.
        mass: unit.Quantity
            Mass of the particles (in units of mass).
        k: unit.Quantity
            Parameter of the external potential (units of energy/length**2) corresponding to the
            stiffness of the harmonic well.

        Examples
        --------

        >>> from uibcdf_test_systems import HarmonicWell
        >>> from openmm import unit
        >>> harmonic_well = HarmonicWell(n_particles=1, mass=32*unit.amu, k=10.0*unit.kilocalorie/(unit.mole*unit.nanometers**2))

        Notes
        -----

        See `corresponding documentation in the user guide regarding this class
        <../../systems/harmonic_well_potential.html>`_.

        """

        super().__init__()

        # Parameters

        self.parameters['n_particles']=n_particles
        self.parameters['mass']=mass
        self.parameters['k']=k

        # OpenMM topology

        self.topology = app.Topology()

        try:
            dummy_element = app.element.get_by_symbol('DUM')
        except:
            dummy_element = app.Element(0, 'DUM', 'DUM', 0.0 * unit.amu)

        chain = self.topology.addChain('A')
        for ii in range(n_particles):
            residue = self.topology.addResidue('DUM', chain)
            atom = self.topology.addAtom(name='DUM', element= dummy_element, residue=residue)

        # OpenMM system

        self.system = mm.System()
        for ii in range(n_particles):
            self.system.addParticle(mass)

        A = 0.5*k # stiffness of the armonic potential for coordinates X, Y and Z

        force = mm.CustomExternalForce('A*(x^2+y^2+z^2)')
        force.addGlobalParameter('A', A)
        for ii in range(n_particles):
            force.addParticle(ii, [])
        _ = self.system.addForce(force)

        # Coordinates

        if coordinates is None:
            coordinates = np.zeros([self.parameters['n_particles'], 3], np.float32) * unit.nanometers

        self.set_coordinates(coordinates)

        # Potential expresion

        x, y, z, k = symbols('x y z k')
        self.potential_expression = (1.0/2.0)*k*(x**2 + y**2 + z**2)
        del(x, y, z, k)

    def evaluate_potential(self, coordinates):

        """Potential evaluation

        The potential energy is evaluated at the position/s specified by the input argument
        `coordinates`.

        Parameters
        ----------

        coordinates: unit.Quantity
            Spatial coordinates of the point or points where the potential energy is evaluated. A
            list, tuple or numpy.ndarray can be used of shape (3) or (n_points,3) with length
            units.

        Returns
        -------

        unit.Quantity
            Value of the energy at the point or points given by the input argment `coordinates`.
            The value of the unit.Quantity will be a single float number or a numpy.ndarray of
            float numbers depending on the shape of `coordinates`.

        Examples
        --------

        >>> from uibcdf_test_systems import HarmonicWell
        >>> from openmm import unit
        >>> harmonic_well = HarmonicWell(n_particles=1, mass=32*unit.amu, k=10.0*unit.kilocalorie/(unit.mole*unit.nanometers**2))
        >>> harmonic_well.potential([-1.5, 0.0, 0.0] * unit.nanometers)
        Quantity(value=11.25, unit=kilocalorie/mole)

        Notes
        -----

        See `corresponding documentation in the user guide regarding this class
        <../../systems/harmonic_well_potential.html>`_.

        """

        coordinates._value = np.array(coordinates._value)
        k = self.parameters['k']

        if len(coordinates._value.shape)==1 and coordinates._value.shape[0]==3:

            x  = coordinates[0]
            y  = coordinates[1]
            z  = coordinates[2]

            return 0.5*k*(x**2 + y**2 + z**2)

        elif len(coordinates._value.shape)==2 and coordinates._value.shape[1]==3:

            x  = coordinates[:,0]
            y  = coordinates[:,1]
            z  = coordinates[:,2]

            return 0.5*k*(x**2 + y**2 + z**2)

        else:

            raise ValueError('The input argument coordinates needs a specific shape.')

    def get_oscillations_time_period(self):

        mass = self.parameters['mass']
        k = self.parameters['k']

        return 2*np.pi*np.sqrt(mass/k)

    def get_standard_deviation(self, temperature=300*unit.kelvin):

        k = self.parameters['k']
        kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
        return np.sqrt(kB*temperature/k)

