import molsysmt as msm
import simtk.unit as unit

sequence = 'aminoacids3:AceAlaNme'

# vacuum.pdb
output='vacuum.pdb'
print("{}... ".format(output), end =" ")
tmp_item = msm.build_peptide(sequence, forcefield='AMBER96', to_form='openmm.Modeller', verbose=False)
msm.energy_minimization(tmp_item, forcefield='AMBER96', to_form=output, verbose=False)
print("done")

