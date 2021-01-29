import molsysmt as msm
import simtk.unit as unit

sequence = 'TyrGlyGlyPheMet'

# vacuum.pdb
output='vacuum.pdb'
print("{}... ".format(output), end =" ")
tmp_item = msm.build_peptide(sequence, forcefield='AMBER14', to_form='openmm.Modeller', verbose=False)
msm.terminal_capping(tmp_item)
msm.energy_minimization(tmp_item, forcefield='AMBER14', to_form=output, verbose=False)
print("done")

# octahedral_14.pdb
output='octahedral_14.pdb'
print("{}... ".format(output, end =" "))
msm.solvate('vacuum.pdb', box_geometry="truncated_octahedral", clearance=14.0*unit.angstroms, water_model='TIP3P',
anion='Cl-', num_anions="neutralize", cation='Na+', num_cations="neutralize",
 18              ionic_strength= 0.0*unit.molar, forcefield='AMBER14', engine="LEaP",
 19              to_form= None, logfile=False, verbose=False):
)
msm.build_peptide(sequence, forcefield='AMBER14', water_model='TIP3P', box_geometry='truncated_octahedral',
                  clearance=14*unit.angstroms, to_form=output, verbose=False)
print("done")

