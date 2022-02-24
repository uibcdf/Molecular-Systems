import openmm.unit as unit

class AlanineTetrapeptideVacuum():

    """Alanine tetrapeptide in vacuum

    Alanine tetrapeptide same as OpenMMTools


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

        >>> from uibcdf_test_systems.systems import AlanineTetrapeptideVacuum
        >>> from openmm import unit
        >>> dialanine = AlanineTetrapeptideVacuum(forcefield='xxx')

        Notes
        -----

        See `corresponding documentation in the user guide regarding this class
        <../../systems/tetraalanine.html>`_.

        """

        # Parameters

        self.parameters={}
        self.parameters['forcefield']=forcefield
        self.parameters['constraints']=constraints
        self.parameters['hydrogen_mass']=hydrogen_mass

        # openmm.Modeller

        from molsysmt import convert
        from .files import alanine_tetrapeptide
        file_pdb = alanine_tetrapeptide['vacuum.pdb']
        item_modeller = convert(file_pdb, to_form='openmm.Modeller')

        # OpenMM topology

        self.topology = item_modeller.topology

        # Coordinates

        self.coordinates = item_modeller.positions

        # OpenMM system

        self.system = convert(item_modeller, to_form='openmm.System', forcefield=forcefield,
                              constraints=constraints, hydrogen_mass=hydrogen_mass,
                              implicit_solvent=None, non_bonded_method='no_cutoff')

class AlanineTetrapeptideImplicitSolvent():

    """Alanine tetrapeptide in implicit solvent

    Alanine tetrapeptide same as OpenMMTools


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

        """Creating a new instance of AlanineTetrapeptideImplicitSolvent

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

        >>> from uibcdf_test_systems.systems import AlanineTetrapeptideImplicitSolvent
        >>> from openmm import unit
        >>> dialanine = AlanineTetrapeptideImplicitSolvent(forcefield='xxx')

        Notes
        -----

        See `corresponding documentation in the user guide regarding this class
        <../../systems/tetraalanine.html>`_.

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
        from .files import alanine_tetrapeptide
        file_pdb = alanine_tetrapeptide['vacuum.pdb']
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

class AlanineTetrapeptideExplicitSolvent():

    """Alanine tetrapeptide in explicit solvent

    Alanine tetrapeptide same as OpenMMTools


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

        """Creating a new instance of AlanineTetrapeptideExplicitSolvent

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

        >>> from uibcdf_test_systems.systems import AlanineTetrapeptideExplicitSolvent
        >>> from openmm import unit
        >>> dialanine = AlanineTetrapeptideExplicitSolvent()

        Notes
        -----

        See `corresponding documentation in the user guide regarding this class
        <../../systems/tetraalanine.html>`_.

        """

        # openmm.Modeller

        from molsysmt import convert
        from .files import alanine_tetrapeptide
        file_pdb = alanine_tetrapeptide['octahedral_14.pdb']
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

