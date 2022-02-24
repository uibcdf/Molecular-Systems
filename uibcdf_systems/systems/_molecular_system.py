import numpy as np
import openmm.unit as unit

class MolecularSystem():

    def __init__(self):

        self.parameters={}
        self.topology=None
        self.system=None
        self.coordinates=None
        self.velocities=None
        self.box=None

    def set_coordinates(self, coordinates):

        coordinates = coordinates.in_units_of(unit.nanometers)
        self.coordinates = np.array(coordinates._value)*unit.nanometers

        pass

    def set_velocities(self, coordinates):

        velocities = velocities.in_units_of(unit.nanometers/unit.picoseconds)
        self.velocities = np.array(velocities._value)*unit.nanometers/unit.picoseconds

        pass

    def set_box(self, box):

        box = box.in_units_of(unit.nanometers)
        self.box = np.array(box._value)*unit.nanometers

        if self.topology is not None:
            self.topology.setPeriodicBoxVectors(box)

        if self.system is not None:
            self.system.setDefaultPeriodicBoxVectors(box[0], box[1], box[2])

        pass

