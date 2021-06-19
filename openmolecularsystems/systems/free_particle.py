from .openmolecularsystem import OpenMolecularSystem
import simtk.unit as unit
import simtk.openmm as mm
import simtk.openmm.app as app
import numpy as np

class FreeParticle(OpenMolecularSystem):

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

    def __init__(self, n_particles=1, mass=16*unit.amus, coordinates=None, box=None):

        """Creating a new instance of FreeParticle

        A new test system is returned with the openmm system of particles behaving freely -no
        potential external potential-.

        Parameters
        ----------

        n_particles: int
            Number of particles in the system
        mass: unit.Quantity
            Mass of the particles (in unites of mass) # mass of CH4

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

        super().__init__()

        # Parameters

        self.parameters['n_particles']=n_particles
        self.parameters['mass']=mass

        # OpenMM topology

        self.topology = None

        # OpenMM system

        self.system = mm.System()

        for ii in range(n_particles):
            self.system.addParticle(mass)

        # Coordinates

        if coordinates is None:
            coordinates = np.zeros([self.parameters['n_particles'], 3], np.float32) * unit.nanometers

        self.set_coordinates(coordinates)


        # Box

        if box is not None:
            self.set_box(box)

