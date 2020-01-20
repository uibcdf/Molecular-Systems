

class FreeParticle():

    """Free particles

    Free particles test system

    Attributes
    ----------
    system : openmm.system
        Openmm system

    """

    system = None

    def __init__(self, n_particles, mass):

       """Creating a new instance of FreeParticle

       A new test system is returned with the openmm system of particles behaving freely -no
       potential external potential-.

       Parameters
       ----------

       n_particles: int
           Number of particles in the system
       mass: unit.Quantity
           Mass of the particles

       Examples
       --------

       >>> from uibcdf_test_systems import FreeParticle
       >>> from simtk import unit
       >>> free_particles = FreeParticle(n_particles=10, mass=32*unit.amu)

       See Also
       --------

       Notes
       -----

       See `the free particle documentation in the user guide section
       <../../systems/free_particle.html>`_.

       """

       # OpenMM system

       import simtk.openmm as mm
       import simtk.unit as unit
       import simtk.openmm.app as app

       self.system = mm.System()

       for ii in range(n_particles):
           self.system.addParticle(mass)

