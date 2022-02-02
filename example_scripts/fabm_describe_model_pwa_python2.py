#!/usr/bin/env python

import sys

if len(sys.argv)!=2:
   print 'This script takes one argument: the path to a YAML file with FABM settings (typically fabm.yaml).'
   sys.exit(2)
yamlfile = sys.argv[1]

try:
   import pyfabm
except ImportError:
   print 'Unable to load pyfabm. Please build and install FABM with FABMHOST=python.'
   sys.exit(1)

# Create model object from YAML file.
model = pyfabm.Model(yamlfile)

# List pelagic state variables
print ''
print 'Pelagic state variables:'
for variable in model.bulk_state_variables:
   print variable.output_name,variable.units,variable.long_path

# List surface state variables
print ''
print 'Surface state variables:'
for variable in model.surface_state_variables:
   print variable.output_name,variable.units,variable.long_path

# List bottom state variables
print ''
print 'Bottom state variables:'
for variable in model.bottom_state_variables:
   print variable.output_name,variable.units,variable.long_path

# List state variable dependencies
print ''
print 'State variable dependencies (forcings):'
for variable in model.dependencies:
   print variable.output_name,variable.units,variable.long_path

# List diagnostics (can be too big for models like ERSEM)
print ''
print 'Bulk model diagnostics:'
#for variable in model.diagnostic_variables:
#   print variable.output_name,variable.units,variable.long_path
for a, variable in enumerate(model.bulk_diagnostic_variables,1):
   print a,variable.output_name,variable.units,variable.long_path

# List horizontal diagnostics
print ''
print 'Horizontal model diagnostics:'
#for variable in model.diagnostic_variables:
#   print variable.output_name,variable.units,variable.long_path
for a, variable in enumerate(model.horizontal_diagnostic_variables,1):
   print a,variable.output_name,variable.units,variable.long_path

# List conserved quantities
print ''
print 'Conserved quantities'
for variable in model.conserved_quantities:
   print variable.name,variable.units

