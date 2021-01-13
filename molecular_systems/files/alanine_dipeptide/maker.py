import molsysmt as msm
import simtk.unit as unit

sequence = 'aminoacids3:AceAlaNme'

# vacuum.pdb
output='vacuum.pdb'
print("{}... ".format(output), end =" ")
msm.build_peptide(sequence, forcefield='AMBER96', to_form=output, verbose=False)
print("done")

# amber96_gbsaobc.prmtop and amber96_gbsaobc.inpcrd
output=['amber96_gbsaobc.prmtop', 'amber96_gbsaobc.inpcrd']
print("{} and {}... ".format(output[0], output[1]), end =" ")
msm.build_peptide(sequence, forcefield='AMBER96',
                  implicit_solvent='OBC1',
                  to_form=output, verbose=False)
print("done")

# amber14_gbsaobc.prmtop and amber14_gbsaobc.inpcrd
output=['amber14_gbsaobc.prmtop', 'amber14_gbsaobc.inpcrd']
print("{} and {}... ".format(output[0], output[1]), end =" ")
msm.build_peptide(sequence, forcefield='AMBER14',
                  implicit_solvent='OBC1',
                  to_form=output, verbose=False)
print("done")

# amber96_tip3p_octahedral_14.prmtop and amber96_tip3p_octahedral_14.inpcrd
output=['amber96_tip3p_octahedral_14.prmtop', 'amber96_tip3p_octahedral_14.inpcrd']
print("{} and {}... ".format(output[0], output[1]), end =" ")
msm.build_peptide(sequence, forcefield='AMBER96',
                  water_model='TIP3P', box_geometry='truncated_octahedral',
                  clearance=14*unit.angstroms,
                  to_form=output, verbose=False)
print("done")


# amber14_tip3p_octahedral_14.prmtop and amber14_tip3p_octahedral_14.inpcrd
output=['amber14_tip3p_octahedral_14.prmtop', 'amber14_tip3p_octahedral_14.inpcrd']
print("{} and {}... ".format(output[0], output[1]), end =" ")
msm.build_peptide(sequence, forcefield='AMBER14',
                  water_model='TIP3P', box_geometry='truncated_octahedral',
                  clearance=14*unit.angstroms,
                  to_form=output, verbose=False)
print("done")

# octahedral_14.pdb
output='octahedral_14.pdb'
print("{}... ".format(output, end =" "))
msm.build_peptide(sequence, forcefield='AMBER14',
                  water_model='TIP3P', box_geometry='truncated_octahedral',
                  clearance=14*unit.angstroms,
                  to_form=output, verbose=False)
print("done")

