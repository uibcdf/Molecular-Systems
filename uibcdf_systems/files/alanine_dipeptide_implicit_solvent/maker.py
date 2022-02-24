import molsysmt as msm
import openmm.unit as unit

sequence = 'aminoacids3:AceAlaNme'

# amber96_gbsaobc.prmtop and amber96_gbsaobc.inpcrd
output=['amber96_gbsaobc.prmtop', 'amber96_gbsaobc.inpcrd']
print("{} and {}... ".format(output[0], output[1]), end =" ")
msm.build_peptide(sequence, forcefield='AMBER96', implicit_solvent='OBC1', to_form=output, verbose=False)
msm.energy_minimization(output, forcefield='AMBER96', implicit_solvent='OBC1', to_form=output[1], verbose=False)
print("done")

# amber14_gbsaobc.prmtop and amber14_gbsaobc.inpcrd
output=['amber14_gbsaobc.prmtop', 'amber14_gbsaobc.inpcrd']
print("{} and {}... ".format(output[0], output[1]), end =" ")
msm.build_peptide(sequence, forcefield='AMBER14', implicit_solvent='OBC1', to_form=output, verbose=False)
msm.energy_minimization(output, forcefield='AMBER14', implicit_solvent='OBC1', to_form=output[1], verbose=False)
print("done")

