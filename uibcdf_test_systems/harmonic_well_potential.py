
class HarmonicWell():

    """Particles in an harmonic well potential

    Test system with particles in an harmonic well potential.

    Attributes
    ----------
    system : openmm.system
        Openmm system
    potential : sympy.function
        External potential expression as a sympy function.

    """


    system = None
    potential = None

    def __init__(self, n_particles, mass, k):

       """Creating a new instance of HarmonicWell

       A new test system is returned with the openmm system of particles in an external harmonic
       well potential:

           $0.5*k*x**2$

       Parameters
       ----------

       n_particles: int
           Number of particles in the system
       mass: unit.Quantity
           Mass of the particles
       k: unit.Quantity
           Parameter of the external potential (units of energy/length**2) corresponding to the
           stiffness of the harmonic well.

       Examples
       --------

       >>> from uibcdf_test_systems import HarmonicWell
       >>> from simtk import unit
       >>> harmonic_well_potential = HarmonicWell(n_particles=1, mass=32*unit.amu, k=10.0*unit.kilocalorie/(unit.mole*unit.nanometers**2))

       See Also
       --------

       """

       ### mass -> unit.amu
       ### m -> unit.kilocalories_per_mole

       # OpenMM system

       import simtk.openmm as mm
       import simtk.unit as unit
       import simtk.openmm.app as app

       self.system = mm.System()

       for ii in range(n_particles):
           self.system.addParticle(mass)

       A = 0.5*k # stiffness of the armonic potential for coordinates X, Y and Z

       force = mm.CustomExternalForce('A*(x^2+y^2+z^2)')
       force.addGlobalParameter('A', A)
       for ii in range(n_particles):
           force.addParticle(ii, [])
       self.system.addForce(force)

       # Potential expresion

       import sympy as sy

       x, y, z, k = sy.symbols('x y z k')
       self.potential = k*(x**2 + y**2 + z**2)
       del(x, y, z, k)


