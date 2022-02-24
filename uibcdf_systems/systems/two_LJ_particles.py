from ._molecular_system import MolecularSystem
import openmm.unit as unit
import openmm as mm
import openmm.app as app
import numpy as np
from sympy import symbols

atoms_LJ={
        'Ar':
             {'mass':39.948*unit.amu,
              'sigma':3.404*unit.angstroms,
              'epsilon':0.238*unit.kilocalories_per_mole
              },
        'Xe':
             {'mass':131.293*unit.amu,
              'sigma':3.961*unit.angstroms,
              'epsilon':0.459*unit.kilocalories_per_mole
              },
        'Dummy_A':
             {'mass':40.0*unit.amu,
              'sigma':2.0*unit.angstroms,
              'epsilon':2.0*unit.kilocalories_per_mole
              },
        'Dummy_B':
             {'mass':120.0*unit.amu,
              'sigma':4.0*unit.angstroms,
              'epsilon':4.5*unit.kilocalories_per_mole
              },
           }

class TwoLJParticles(MolecularSystem):

    def __init__(self, atom_1='Ar', atom_2='Xe', mass_1=None, sigma_1=None, epsilon_1=None,
                 mass_2=None, sigma_2=None, epsilon_2=None,
                 cutoff_distance = None, switching_distance = None, box=None,
                 coordinates=None):

        super().__init__()

        # Parameters

        if (mass_1 is not None) or (sigma_1 is not None) or (epsilon_1 is not None):
            if mass_1 is None:
                raise ValueError('A value for the input argument "mass_1" is needed.')
            if sigma_1 is None:
                raise ValueError('A value for the input argument "sigma_1" is needed.')
            if epsilon_1 is None:
                raise ValueError('A value for the input argument "epsilon_1" is needed.')
        elif atom_1 is not None:
            mass_1 = atoms_LJ[atom_1]['mass']
            sigma_1 = atoms_LJ[atom_1]['sigma']
            epsilon_1 = atoms_LJ[atom_1]['epsilon']

        if (mass_2 is not None) or (sigma_2 is not None) or (epsilon_2 is not None):
            if mass_2 is None:
                raise ValueError('A value for the input argument "mass_2" is needed.')
            if sigma_2 is None:
                raise ValueError('A value for the input argument "sigma_2" is needed.')
            if epsilon_2 is None:
                raise ValueError('A value for the input argument "epsilon_2" is needed.')
        elif atom_2 is not None:
            mass_2 = atoms_LJ[atom_2]['mass']
            sigma_2 = atoms_LJ[atom_2]['sigma']
            epsilon_2 = atoms_LJ[atom_2]['epsilon']

        self.parameters['mass_1']=mass_1
        self.parameters['sigma_1']=sigma_1
        self.parameters['epsilon_1']=epsilon_1

        self.parameters['mass_2']=mass_2
        self.parameters['sigma_2']=sigma_2
        self.parameters['epsilon_2']=epsilon_2

        if box is not None:

            reduced_sigma = self.get_reduced_sigma()

            if cutoff_distance is None:
                cutoff_distance = 4.0*reduced_sigma

            if switching_distance is None:
                switching_distance = 3.0*reduced_sigma

        self.parameters['box'] = box
        self.parameters['cutoff_distance'] = cutoff_distance
        self.parameters['switching_distance'] = switching_distance

        # OpenMM topology

        self.topology = app.Topology()

        try:
            dummy_element = app.element.get_by_symbol('DUM')
        except:
            dummy_element = app.Element(0, 'DUM', 'DUM', 0.0 * unit.amu)

        chain = self.topology.addChain('A')
        residue = self.topology.addResidue('DUM_0', chain)
        atom = self.topology.addAtom(name='DUM_0', element= dummy_element, residue=residue)
        residue = self.topology.addResidue('DUM_1', chain)
        atom = self.topology.addAtom(name='DUM_1', element= dummy_element, residue=residue)

        # OpenMM system

        self.system = mm.System()

        non_bonded_force = mm.NonbondedForce()

        if box is not None:
            non_bonded_force.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
            non_bonded_force.setUseSwitchingFunction(True)
            non_bonded_force.setCutoffDistance(cutoff_distance)
            non_bonded_force.setSwitchingDistance(switching_distance)
        else:
            non_bonded_force.setNonbondedMethod(mm.NonbondedForce.NoCutoff)

        self.system.addParticle(mass_1)
        charge_1 = 0.0 * unit.elementary_charge
        non_bonded_force.addParticle(charge_1, sigma_1, epsilon_1)

        self.system.addParticle(mass_2)
        charge_2 = 0.0 * unit.elementary_charge
        non_bonded_force.addParticle(charge_2, sigma_2, epsilon_2)

        _ = self.system.addForce(non_bonded_force)

        # Coordinates

        if coordinates is not None:
            self.set_coordinates(coordinates)

        # Box

        if box is not None:
            self.set_box(box)

        # Potential expresion

        d, eps_r, sigma_r = symbols('d eps_r sigma_r')
        self.potential_expression = 4.0*eps_r * ((sigma_r/d)**12 - (sigma_r/d)**6)
        del(d, eps_r, sigma_r)

    def get_reduced_sigma(self):

        return 0.5*(self.parameters['sigma_1']+self.parameters['sigma_2'])

    def get_reduced_epsilon(self):

        return np.sqrt(self.parameters['epsilon_1']*self.parameters['epsilon_2'])

    def get_reduced_mass(self):
        return (self.parameters['mass_1']*self.parameters['mass_2'])/(self.parameters['mass_1']+self.parameters['mass_2'])

    def evaluate_potential(self, coordinates):

        coordinates = coordinates.in_units_of(unit.angstroms)
        reduced_epsilon = self.get_reduced_epsilon()
        reduced_epsilon = reduced_epsilon.in_units_of(unit.kilocalories_per_mole)
        reduced_sigma = self.get_reduced_sigma()
        reduced_sigma = reduced_sigma.in_units_of(unit.angstroms)

        dx  = coordinates[0,0]-coordinates[0,0]
        dy  = coordinates[0,1]-coordinates[0,1]
        dz  = coordinates[0,2]-coordinates[0,2]
        d = np.sqrt(dx**2+dy**2+dz**2)
        c = reduced_sigma/d

        return 4.0*reduced_epsilon*(c**12-c**6)

    def get_coordinates_minimum(self):

        reduced_sigma = self.get_reduced_sigma()

        return 2.0**(1.0/6.0) * reduced_sigma


    def get_small_oscillations_time_period_around_minimum(self):

        reduced_mass = self.get_reduced_mass()
        reduced_sigma = self.get_reduced_sigma()
        reduced_epsilon = self.get_reduced_epsilon()
        k = 36.0*2.0**(2.0/3.0)*reduced_epsilon/reduced_sigma**2
        tau = 2.0*np.pi*np.sqrt(reduced_mass/k)

        return tau

    def get_harmonic_standard_deviation_around_minimum(self, temperature=300*unit.kelvin):

        reduced_mass = self.get_reduced_mass()
        reduced_sigma = self.get_reduced_sigma()
        reduced_epsilon = self.get_reduced_epsilon()
        k = 36.0*2.0**(2.0/3.0)*reduced_epsilon/reduced_sigma**2
        tau = 2.0*np.pi*np.sqrt(reduced_mass/k)

        kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
        return np.sqrt(kB*temperature/k)

