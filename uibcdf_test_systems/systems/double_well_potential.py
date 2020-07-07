from simtk.unit import amu, kilocalories_per_mole, nanometers, picoseconds
from numpy import zeros

class DoubleWell():

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

    n_particles = None
    mass = None
    potential_expression = None
    potential_parameters = None

    forcefield = None
    system_parameters = None

    system = None
    coordinates = None
    topology = None


    def __init__(self, n_particles=1, mass=64*amu, Eo=4.0*kilocalories_per_mole, a=1.0*nanometers,
                 b=0.0*kilocalories_per_mole, coordinates= zeros([1,3],dtype=float)*nanometers):

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

        # OpenMM system

        import simtk.openmm as mm
        import simtk.unit as unit
        import simtk.openmm.app as app

        self.system = mm.System()
        self.n_particles = n_particles
        self.mass = mass
        self.coordinates = coordinates
        self.topology = None

        for ii in range(n_particles):
            self.system.addParticle(mass)

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
        self.system.addForce(force)

        # Potential expresion and constants

        from sympy import symbols

        self.potential_parameters ={'Eo':Eo, 'a':a, 'b':b}

        x, y, z, Eo, a, b = symbols('x y z Eo a b')
        self.potential_expression = Eo*((x/a)**4-2.0*(x/a)**2)-(b/a)*x + 0.5 *(8.0*Eo/a**2)*(y**2 + z**2)
        del(x, y, z, Eo, a, b)

    def potential_energy(self, coordinates=None):

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

        from numpy import array

        Eo = self.potential_parameters['Eo']
        a = self.potential_parameters['a']
        b = self.potential_parameters['b']

        if coordinates is None:
            coordinates = self.coordinates
        else:
            coordinates._value = array(coordinates._value)

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

    def coordinates_minima(self):

        import sympy as sy

        Eo = self.potential_parameters['Eo']
        a = self.potential_parameters['a']
        b = self.potential_parameters['b']

        x, y, z = sy.symbols('x y z')
        xu = x*nanometers
        yu = y*nanometers
        zu = z*nanometers
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
                root_3d=zeros([3],dtype=float)*nanometers
                root_3d[0]=root*nanometers
                roots.append(root_3d)

        del(x, y, z)

        return roots

    def coordinates_maxima(self):

        import sympy as sy

        Eo = self.potential_parameters['Eo']
        a = self.potential_parameters['a']
        b = self.potential_parameters['b']

        x, y, z = sy.symbols('x y z')
        xu = x*nanometers
        yu = y*nanometers
        zu = z*nanometers

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
                root_3d=zeros([3],dtype=float)*nanometers
                root_3d[0]=root*nanometers
                roots.append(root_3d)

        del(x, y, z)

        return roots

    def armonic_oscillation_periods(self):

        import sympy as sy
        import numpy as np

        Eo = self.potential_parameters['Eo']
        a = self.potential_parameters['a']
        b = self.potential_parameters['b']

        x, y, z = sy.symbols('x y z')
        xu = x*nanometers
        yu = y*nanometers
        zu = z*nanometers

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
                root_y=root*nanometers
                T_y = 2*np.pi*np.sqrt(self.mass/(effective_k * kilocalories_per_mole/nanometers**2))

        g=sy.diff(potential_x,x)
        gg=sy.diff(potential_x,x,x)
        roots_diff=sy.roots(g,x)

        roots=[]
        Ts=[]
        for root in roots_diff.keys():
            effective_k=gg.subs(x,root)
            if effective_k>0:
                root_3d=zeros([3],dtype=float)*nanometers
                root_3d[0]=root*nanometers
                root_3d[1]=root_y
                root_3d[2]=root_y
                roots.append(root_3d)
                T_3d=zeros([3],dtype=float)*picoseconds
                T_3d[0] = 2*np.pi*np.sqrt(self.mass/(effective_k * kilocalories_per_mole/nanometers**2))
                T_3d[1] = T_y
                T_3d[2] = T_y
                Ts.append(T_3d)

        del(x, y, z)

        return roots, Ts

