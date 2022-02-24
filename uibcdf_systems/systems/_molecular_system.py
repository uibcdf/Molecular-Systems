from uibcdf_systems import puw
from uibcdf_systems.tools import energy_minimization as _energy_minimization
import openmm as mm
import openmm.unit as unit
import numpy as np

class MolecularSystem():

    def __init__(self):

        self.parameters={}
        self.topology=None
        self.system=None
        self.coordinates=None
        self.velocities=None
        self.box=None

    def set_coordinates(self, coordinates):

        self.coordinates = puw.convert(coordinates, to_unit='nm')

        pass

    def set_velocities(self, velocities):

        self.velocities = puw.convert(velocities, to_unit='nm/ps')

        pass

    def set_box(self, box):

        self.box = puw.convert(box, unit='nm')

        if self.topology is not None:
            self.topology.setPeriodicBoxVectors(box)

        if self.system is not None:
            self.system.setDefaultPeriodicBoxVectors(box[0], box[1], box[2])

        pass

    def energy_minimization(self, platform_name='CUDA', verbose=False):

        _energy_minimization(self, plaform_name=platform_name, verbose=verbose)

        pass

