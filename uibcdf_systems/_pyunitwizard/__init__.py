import pyunitwizard as puw
puw.configure.load_library(['openmm.unit', 'pint'])
puw.configure.set_default_form('openmm.unit')
puw.configure.set_default_parser('pint')
puw.configure.set_standard_units(['nm', 'ps', 'K', 'mole', 'amu', 'e',
                                 'kJ/mol', 'kJ/(mol*nm**2)', 'N', 'degrees'])

