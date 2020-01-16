
class DoubleWell():

    system = None
    potential_expression = None
    potential_parameters = None

    def __init__(self, n_particles, mass, Eo, a, b):

        ### mass -> unit.amu
        ### Eo -> unit.kilocalories_per_mole
        ### a -> unit.nanometers
        ### b -> unit.kilocalories_per_mole

        # OpenMM system

        import simtk.openmm as mm
        import simtk.unit as unit
        import simtk.openmm.app as app

        self.system = mm.System()

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
        self.potential._expression = Eo*((x/a)**4-2.0*(x/a)**2)-(b/a)*x + 0.5 *(8.0*Eo/a**2)*(y**2 + z**2)
        del(x, y, z, Eo, a, b)

    def potential(self):
        pass

