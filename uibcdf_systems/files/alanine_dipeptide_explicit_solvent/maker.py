import molsysmt as msm
import openmm.unit as unit

sequence = 'aminoacids3:AceAlaNme'

# amber96_tip3p_octahedral_14.prmtop and amber96_tip3p_octahedral_14.inpcrd
output=['amber96_tip3p_octahedral_14.prmtop', 'amber96_tip3p_octahedral_14.inpcrd']
print("{} and {}... ".format(output[0], output[1]), end =" ")
msm.build_peptide(sequence, forcefield='AMBER96', water_model='TIP3P', box_geometry='truncated_octahedral',
                  clearance=14*unit.angstroms, to_form=output, verbose=False)
msm.energy_minimization(output, forcefield='AMBER96', water_model='TIP3P', to_form=output[1], verbose=False)
print("done")

# amber14_tip3p_octahedral_14.prmtop and amber14_tip3p_octahedral_14.inpcrd
output=['amber14_tip3p_octahedral_14.prmtop', 'amber14_tip3p_octahedral_14.inpcrd']
print("{} and {}... ".format(output[0], output[1]), end =" ")
msm.build_peptide(sequence, forcefield='AMBER14', water_model='TIP3P', box_geometry='truncated_octahedral',
                  clearance=14*unit.angstroms, to_form=output, verbose=False)
msm.energy_minimization(output, forcefield='AMBER14', water_model='TIP3P', to_form=output[1], verbose=False)
print("done")

# octahedral_14.pdb
output='octahedral_14.pdb'
print("{}... ".format(output, end =" "))
tmp_item = msm.build_peptide(sequence, forcefield='AMBER14', water_model='TIP3P', box_geometry='truncated_octahedral',
                  clearance=14*unit.angstroms, to_form='openmm.Modeller', verbose=False)
msm.energy_minimization(output, forcefield='AMBER14', water_model='TIP3P', to_form=output, verbose=False)
print("done")

