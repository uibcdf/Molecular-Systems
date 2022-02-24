import openmm.unit as unit

class AlanineDipeptideVacuum():

    """Alanine dipeptide in vacuum

    Alanine dipeptide same as OpenMMTools


    Attributes
    ----------
    system
        Xxx xxx

    Methods
    -------

    """

    def __init__(self, forcefield='AMBER14', constraints='h_bonds', hydrogen_mass=None):

        """XXX

        A new test system is returned with the openmm system of particles in an external double
        well potential.

        Parameters
        ----------

        forcefield: xxx
            XXX
        constraints: app.Hbonds
            XXX
        hydrogen_mass: None
            XXX

        Examples
        --------

        >>> from uibcdf_test_systems.systems import AlanineDipeptideVacuum
        >>> from openmm import unit
        >>> dialanine = AlanineDipeptideVacuum(forcefield='xxx')

        Notes
        -----

        See `corresponding documentation in the user guide regarding this class
        <../../systems/dialanine.html>`_.

        """

        # Parameters

        self.parameters={}
        self.parameters['forcefield']=forcefield
        self.parameters['constraints']=constraints
        self.parameters['hydrogen_mass']=hydrogen_mass

        # openmm.Modeller

        from molsysmt import convert
        from .files import alanine_dipeptide
        file_pdb = alanine_dipeptide['vacuum.pdb']
        item_modeller = convert(file_pdb, to_form='openmm.Modeller')

        # OpenMM topology

        self.topology = item_modeller.topology

        # Coordinates

        self.coordinates = item_modeller.positions

        # OpenMM system

        self.system = convert(item_modeller, to_form='openmm.System', forcefield=forcefield,
                              constraints=constraints, hydrogen_mass=hydrogen_mass,
                              implicit_solvent=None, non_bonded_method='no_cutoff')

class AlanineDipeptideImplicitSolvent():

    """Alanine dipeptide in implicit solvent

    Alanine dipeptide same as OpenMMTools


    Attributes
    ----------
    system
        Xxx xxx
    positions
        Xxx xxx
    topology
        Xxx xxx

    Methods
    -------

    """

    def __init__(self, forcefield='AMBER14', implicit_solvent='OBC1', solute_dielectric=1.0,
                 solvent_dielectric=78.3, implicit_solvent_kappa=0.0/unit.nanometer,
                 implicit_solvent_salt_conc=0.0*unit.mole/unit.liter,
                 constraints='h_bonds', hydrogen_mass=None):

        """Creating a new instance of AlanineDipeptideImplicitSolvent

        ...

        Parameters
        ----------

        forcefield: xxx
            XXX
        constraints: app.Hbonds
            XXX
        hydrogen_mass: None
            XXX

        Examples
        --------

        >>> from uibcdf_test_systems.systems import AlanineDipeptideImplicitSolvent
        >>> from openmm import unit
        >>> dialanine = AlanineDipeptideImplicitSolvent(forcefield='xxx')

        Notes
        -----

        See `corresponding documentation in the user guide regarding this class
        <../../systems/dialanine.html>`_.

        """

        # Parameters

        self.parameters={}
        self.parameters['forcefield']=forcefield
        self.parameters['implicit_solvent']=implicit_solvent
        self.parameters['solute_dielectric']=solute_dielectric
        self.parameters['solvent_dielectric']=solvent_dielectric
        self.parameters['implicit_solvent_kappa']=implicit_solvent_kappa
        self.parameters['implicit_solvent_salt_conc']=implicit_solvent_salt_conc
        self.parameters['constraints']=constraints
        self.parameters['hydrogen_mass']=hydrogen_mass

        # openmm.Modeller

        from molsysmt import convert
        from .files import alanine_dipeptide
        file_pdb = alanine_dipeptide['vacuum.pdb']
        item_modeller = convert(file_pdb, to_form='openmm.Modeller')

        # OpenMM topology

        self.topology = item_modeller.topology

        # Coordinates

        self.coordinates = item_modeller.positions

        # OpenMM system

        self.system = convert(item_modeller, to_form='openmm.System', forcefield=forcefield,
                              constraints=constraints, hydrogen_mass=hydrogen_mass, non_bonded_method='no_cutoff',
                              implicit_solvent=implicit_solvent, solute_dielectric=solute_dielectric,
                              solvent_dielectric=solvent_dielectric,
                              implicit_solvent_kappa=implicit_solvent_kappa,
                              implicit_solvent_salt_conc=implicit_solvent_salt_conc
                             )

class AlanineDipeptideExplicitSolvent():

    """Alanine dipeptide in explicit solvent

    Alanine dipeptide same as OpenMMTools


    Attributes
    ----------
    system
        Xxx xxx
    positions
        Xxx xxx
    topology
        Xxx xxx

    Methods
    -------

    """

    def __init__(self, forcefield = 'AMBER14', water_model = 'TIP3P', rigid_water = True, constraints = 'h_bonds', hydrogen_mass = None,
                 non_bonded_method = 'PME', non_bonded_cutoff = 10.0 * unit.angstroms, switch_distance = 8.5 * unit.angstroms,
                 use_dispersion_correction = True, ewald_error_tolerance = 1.0e-5):

        """Creating a new instance of AlanineDipeptideExplicitSolvent

        Parameters
        ----------

        forcefield: xxx
            XXX
        constraints: app.Hbonds
            XXX
        hydrogen_mass: None
            XXX

        Examples
        --------

        >>> from uibcdf_test_systems.systems import AlanineDipeptideExplicitSolvent
        >>> from openmm import unit
        >>> dialanine = AlanineDipeptideExplicitSolvent()

        Notes
        -----

        See `corresponding documentation in the user guide regarding this class
        <../../systems/dialanine.html>`_.

        """

        # openmm.Modeller

        from molsysmt import convert
        from .files import alanine_dipeptide
        file_pdb = alanine_dipeptide['octahedral_14.pdb']
        item_modeller = convert(file_pdb, to_form='openmm.Modeller')

        # OpenMM topology

        self.topology = item_modeller.topology

        # Coordinates

        self.coordinates = item_modeller.positions

        # OpenMM system

        self.system = convert(item_modeller, to_form='openmm.System', forcefield=forcefield,
                              water_model=water_model, rigid_water=rigid_water, constraints=constraints,
                              hydrogen_mass=hydrogen_mass, non_bonded_method=non_bonded_method,
                              non_bonded_cutoff=non_bonded_cutoff, switch_distance=switch_distance,
                              use_dispersion_correction=use_dispersion_correction,
                              ewald_error_tolerance=ewald_error_tolerance)

