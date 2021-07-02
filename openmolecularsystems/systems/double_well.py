from .openmolecularsystem import OpenMolecularSystem
import simtk.unit as unit
import numpy as np
import sympy as sy
import simtk.unit as unit
import simtk.openmm as mm
import simtk.openmm.app as app


class DoubleWell(OpenMolecularSystem):

    """Particles in an double well potential

    Test system with particles in a quadratic double well potential.

    .. math::

        Eo\\left[\\left(\\frac{x}{a}\\right)^4-2\\left(\\frac{x}{a}\\right)^2\\right]-\\frac{b}{a}x + \\frac{4Eo}{a^{2}}\\left(y^2 + z^2\\right)

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

    def __init__(self, n_particles=1, mass=64*unit.amu, Eo=4.0*unit.kilocalories_per_mole,
                 a=1.0*unit.nanometers, b=0.0*unit.kilocalories_per_mole,
                 coordinates= None):

        """Creating a new instance of DoubleWell

        A new test system is returned with the openmm system of particles in an external double
        well potential.

        Parameters
        ----------

        n_particles: int
            Number of particles in the system
        mass: unit.Quantity
            Mass of the particles (in units of mass).
        Eo: unit.Quantity
            Parameter of the external potential with units of energy.
        a: unit.Quantity
            Parameter of the external potential with units of length.
        b: unit.Quantity
            Parameter of the external potential with units of energy.

        Examples
        --------

        >>> from uibcdf_test_systems import DoubleWell
        >>> from simtk import unit
        >>> double_well = DoubleWell(n_particles = 1, mass = 64 * unit.amu, Eo=4.0 * unit.kilocalories_per_mole, a=1.0 * unit.nanometers, b=0.0 * unit.kilocalories_per_mole))

        Notes
        -----

        See `corresponding documentation in the user guide regarding this class
        <../../systems/double_well_potential.html>`_.

        """

        super().__init__()

        # Parameters

        self.parameters={}
        self.parameters['n_particles']=n_particles
        self.parameters['mass']=mass
        self.parameters['Eo']=Eo
        self.parameters['a']=a
        self.parameters['b']=b

        # OpenMM topology

        self.topology = app.Topology()

        try:
            dummy_element = app.element.get_by_symbol('DUM')
        except:
            dummy_element = app.Element(0, 'DUM', 'DUM', 0.0 * unit.amu)

        dummy_element.mass._value = mass.value_in_unit(unit.amu)

        chain = self.topology.addChain('A')
        for _ in range(n_particles):
            residue = self.topology.addResidue('DUM', chain)
            atom = self.topology.addAtom(name='DUM', element= dummy_element, residue=residue)

        # OpenMM system

        self.system = mm.System()

        for _ in range(n_particles):
            self.system.addParticle(dummy_element.mass)

        k = 8.0*Eo/(a**2) # stiffness of the armonic potential for coordinates Y and Z

        A = Eo/(a**4)
        B = -2.0*Eo/(a**2)
        C = -b/a
        D = k/2.0

        force = mm.CustomExternalForce('A*x^4+B*x^2+C*x + D*(y^2+z^2)')
        force.addGlobalParameter('A', A)
        force.addGlobalParameter('B', B)
        force.addGlobalParameter('C', C)
        force.addGlobalParameter('D', D)

        for ii in range(n_particles):
            force.addParticle(ii, [])
        _ = self.system.addForce(force)

        # Coordinates

        if coordinates is None:
            coordinates = np.zeros([self.parameters['n_particles'], 3], np.float32) * unit.nanometers

        self.set_coordinates(coordinates)

        # Potential expresion and constants

        x, y, z, Eo, a, b = sy.symbols('x y z Eo a b')
        self.potential_expression = Eo*((x/a)**4-2.0*(x/a)**2)-(b/a)*x + 0.5 *(8.0*Eo/a**2)*(y**2 + z**2)
        del(x, y, z, Eo, a, b)

    def evaluate_potential(self, coordinates=None):

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

        >>> from uibcdf_test_systems import DoubleWell
        >>> from simtk import unit
        >>> double_well = DoubleWell(n_particles = 1, mass = 64 * unit.amu, Eo=4.0 * unit.kilocalories_per_mole, a=1.0 * unit.nanometers, b=0.0 * unit.kilocalories_per_mole))
        >>> double_well.potential([-1.5, 0.0, 0.0] * unit.nanometers)
        Quantity(value=2.25, unit=kilocalorie/mole)

        Notes
        -----

        See `corresponding documentation in the user guide regarding this class
        <../../systems/double_well_potential.html>`_.

        """

        Eo = self.parameters['Eo']
        a = self.parameters['a']
        b = self.parameters['b']

        if coordinates is None:
            coordinates = self.coordinates
        else:
            coordinates._value = np.array(coordinates._value)

        if len(coordinates._value.shape)==1 and coordinates._value.shape[0]==3:

            x  = coordinates[0]
            y  = coordinates[1]
            z  = coordinates[2]

            return Eo*((x/a)**4-2.0*(x/a)**2)-(b/a)*x + 0.5 *(8.0*Eo/a**2)*(y**2 + z**2)

        elif len(coordinates._value.shape)==2 and coordinates._value.shape[1]==3:

            x  = coordinates[:,0]
            y  = coordinates[:,1]
            z  = coordinates[:,2]

            return Eo*((x/a)**4-2.0*(x/a)**2)-(b/a)*x + 0.5 *(8.0*Eo/a**2)*(y**2 + z**2)

        else:

            raise ValueError('The input argument coordinates needs a specific shape.')

    def get_coordinates_minima(self):

        Eo = self.parameters['Eo']
        a = self.parameters['a']
        b = self.parameters['b']

        x, y, z = sy.symbols('x y z')
        xu = x*unit.nanometers
        yu = y*unit.nanometers
        zu = z*unit.nanometers
        potential_x = Eo*((xu/a)**4-2.0*(xu/a)**2)-(b/a)*xu
        potential_y = 0.5 *(8.0*Eo/a**2)*(yu**2)
        potential_z = 0.5 *(8.0*Eo/a**2)*(zu**2)

        g=sy.diff(potential_x,x)
        gg=sy.diff(potential_x,x,x)
        roots_diff=sy.roots(g,x)

        roots=[]
        for root in roots_diff.keys():
            effective_k=gg.subs(x,root)
            if effective_k>0:
                root_3d=np.zeros([3],dtype=float)*unit.nanometers
                root_3d[0]=root*unit.nanometers
                roots.append(root_3d)

        del(x, y, z)

        return roots

    def get_coordinates_maximum(self):

        Eo = self.parameters['Eo']
        a = self.parameters['a']
        b = self.parameters['b']

        x, y, z = sy.symbols('x y z')
        xu = x*unit.nanometers
        yu = y*unit.nanometers
        zu = z*unit.nanometers

        potential_x = Eo*((xu/a)**4-2.0*(xu/a)**2)-(b/a)*xu
        potential_y = 0.5 *(8.0*Eo/a**2)*(yu**2)
        potential_z = 0.5 *(8.0*Eo/a**2)*(zu**2)

        g=sy.diff(potential_x,x)
        gg=sy.diff(potential_x,x,x)
        roots_diff=sy.roots(g,x)

        roots=[]
        for root in roots_diff.keys():
            effective_k=gg.subs(x,root)
            if effective_k<0:
                root_3d=np.zeros([3],dtype=float)*unit.nanometers
                root_3d[0]=root*unit.nanometers
                roots.append(root_3d)

        del(x, y, z)

        return roots

    def get_small_oscillations_time_periods_around_minima(self):

        Eo = self.parameters['Eo']
        a = self.parameters['a']
        b = self.parameters['b']
        mass = self.parameters['mass']

        x, y, z = sy.symbols('x y z')
        xu = x*unit.nanometers
        yu = y*unit.nanometers
        zu = z*unit.nanometers

        potential_x = Eo*((xu/a)**4-2.0*(xu/a)**2)-(b/a)*xu
        potential_y = 0.5 *(8.0*Eo/a**2)*(yu**2)
        potential_z = 0.5 *(8.0*Eo/a**2)*(zu**2)

        g=sy.diff(potential_y,y)
        gg=sy.diff(potential_y,y,y)
        roots_diff=sy.roots(g,y)

        root_y=None
        T_y=None
        for root in roots_diff.keys():
            effective_k=gg.subs(y,root)
            if effective_k>0:
                root_y=root*unit.nanometers
                T_y = 2*np.pi*np.sqrt(mass/(effective_k * unit.kilocalories_per_mole/unit.nanometers**2))

        g=sy.diff(potential_x,x)
        gg=sy.diff(potential_x,x,x)
        roots_diff=sy.roots(g,x)

        roots=[]
        Ts=[]
        for root in roots_diff.keys():
            effective_k=gg.subs(x,root)
            if effective_k>0:
                root_3d=np.zeros([3],dtype=float)*unit.nanometers
                root_3d[0]=root*unit.nanometers
                root_3d[1]=root_y
                root_3d[2]=root_y
                roots.append(root_3d)
                T_3d=np.zeros([3],dtype=float)*unit.picoseconds
                T_3d[0] = 2*np.pi*np.sqrt(mass/(effective_k * unit.kilocalories_per_mole/unit.nanometers**2))
                T_3d[1] = T_y
                T_3d[2] = T_y
                Ts.append(T_3d)

        del(x, y, z)

        return roots, Ts

    def get_harmonic_standard_deviations_around_minima(self, temperature=300.0*unit.kelvin):


        kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA

        Eo = self.parameters['Eo']
        a = self.parameters['a']
        b = self.parameters['b']

        x, y, z = sy.symbols('x y z')
        xu = x*unit.nanometers
        yu = y*unit.nanometers
        zu = z*unit.nanometers

        potential_x = Eo*((xu/a)**4-2.0*(xu/a)**2)-(b/a)*xu
        potential_y = 0.5 *(8.0*Eo/a**2)*(yu**2)
        potential_z = 0.5 *(8.0*Eo/a**2)*(zu**2)

        g=sy.diff(potential_y,y)
        gg=sy.diff(potential_y,y,y)
        roots_diff=sy.roots(g,y)

        root_y=None
        sigma_y=None
        for root in roots_diff.keys():
            effective_k=gg.subs(y,root)
            if effective_k>0:
                root_y=root*unit.nanometers
                sigma_y = np.sqrt(kB*temperature/(effective_k*unit.kilocalories_per_mole/unit.nanometers**2))

        g=sy.diff(potential_x,x)
        gg=sy.diff(potential_x,x,x)
        roots_diff=sy.roots(g,x)

        roots=[]
        sigmas=[]
        for root in roots_diff.keys():
            effective_k=gg.subs(x,root)
            if effective_k>0:
                root_3d=np.zeros([3],dtype=float)*unit.nanometers
                root_3d[0]=root*unit.nanometers
                root_3d[1]=root_y
                root_3d[2]=root_y
                roots.append(root_3d)
                sigma_3d=np.zeros([3],dtype=float)*unit.nanometers
                sigma_3d[0] = np.sqrt(kB*temperature/(effective_k*unit.kilocalories_per_mole/unit.nanometers**2))
                sigma_3d[1] = sigma_y
                sigma_3d[2] = sigma_y
                sigmas.append(sigma_3d)

        del(x, y, z)

        return roots, sigmas

