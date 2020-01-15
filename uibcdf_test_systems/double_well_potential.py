

class DoubleWell():

    system = None
    potential = None

    def __init__(self, n_particles, mass, Eo, c, m):

        ### mass -> unit.amu
        ### Eo -> unit.kilocalories_per_mole
        ### c -> unit.nanometers
        ### m -> unit.kilocalories_per_mole

        # OpenMM system

        import simtk.openmm as mm
        import simtk.unit as unit
        import simtk.openmm.app as app

        self.system = mm.System()

        for ii in range(n_particles):
            self.system.addParticle(mass)

        k = 8.0*Eo/(c**2) # stiffness of the armonic potential for coordinates Y and Z

        A = Eo/(c**4)
        B = -2.0*Eo/(c**2)
        C = -m/c
        D = k/2.0

        force = mm.CustomExternalForce('A*x^4+B*x^2+C*x + D*(y^2+z^2)')
        force.addGlobalParameter('A', A)
        force.addGlobalParameter('B', B)
        force.addGlobalParameter('C', C)
        force.addGlobalParameter('D', D)
        for ii in range(n_particles):
            force.addParticle(ii, [])
        self.system.addForce(force)

        # Potential expresion

        import sympy as sy

        x, y, z, Eo, c, m = sy.symbols('x y z Eo c m')
        self.potential = Eo*((x/c)**4-2.0*(x/c)**2)-(m/c)*x + 0.5 *(8.0*Eo/c**2)*(y**2 + z**2)
        del(x, y, z, Eo, c, m)

