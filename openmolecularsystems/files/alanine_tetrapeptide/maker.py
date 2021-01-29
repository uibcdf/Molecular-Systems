import molsysmt as msm
import simtk.unit as unit

sequence = 'aminoacids3:AceAlaAlaAlaNme'

# vacuum.pdb
output='vacuum.pdb'
print("{}... ".format(output), end =" ")
tmp_item = msm.build_peptide(sequence, forcefield='AMBER14', to_form='openmm.Modeller', verbose=False)
msm.energy_minimization(tmp_item, forcefield='AMBER14', to_form=output, verbose=False)
print("done")

# octahedral_14.pdb
output='octahedral_14.pdb'
print("{}... ".format(output, end =" "))
msm.build_peptide(sequence, forcefield='AMBER14', water_model='TIP3P', box_geometry='truncated_octahedral',
                  clearance=14*unit.angstroms, to_form=output, verbose=False)
print("done")

