import numpy as np
import simtk.unit as unit

class OpenMolecularSystem():

    def __init__(self):

        self.parameters={}
        self.topology=None
        self.system=None
        self.coordinates=None
        self.box=None

    def set_coordinates(self, coordinates):

        coordinates = coordinates.in_units_of(unit.nanometers)
        self.coordinates = np.array(coordinates._value)*unit.nanometers

        pass

    def set_box(self, box):

        box = box.in_units_of(unit.nanometers)
        self.box = np.array(box._value)*unit.nanometers

        pass

