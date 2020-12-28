import simtk.unit as unit
import simtk.openmm as mm
import simtk.openmm.app as app
import numpy as np

class FreeParticle():

    """Free particles

    Free particles test system

    Attributes
    ----------
    parameters :
        Dictionary with parameters to initializate the molecular system
    topology : openmm.topology
        Openmm topology
    coordinates :
        Position of particles
    system : openmm.system
        Openmm system

    """

    parameters = {}
    topology = None
    coordinates = None
    system = None


    def __init__(self, n_particles=1, mass=16*unit.amus, coordinates=[[0.0, 0.0, 0.0]]*unit.nanometers):  # CH4 as input example

       """Creating a new instance of FreeParticle

       A new test system is returned with the openmm system of particles behaving freely -no
       potential external potential-.

       Parameters
       ----------

       n_particles: int
           Number of particles in the system
       mass: unit.Quantity
           Mass of the particles (in unites of mass)

       Examples
       --------

       >>> from uibcdf_test_systems import FreeParticle
       >>> from simtk import unit
       >>> free_particle = FreeParticle(n_particles=1, mass=16*unit.amu)

       Notes
       -----

       See `the free particle documentation in the user guide section
       <../../systems/free_particle.html>`_.

       """

       # Parameters

       self.parameters['n_particles']=n_particles
       self.parameters['mass']=mass

       # Coordinates

       coordinates = coordinates.in_units_of(unit.nanometers)
       self.coordinates = np.array(coordinates._value)*unit.nanometers

       # OpenMM topology

       self.topology = None

       # OpenMM system

       self.system = mm.System()

       for ii in range(n_particles):
           self.system.addParticle(mass)

