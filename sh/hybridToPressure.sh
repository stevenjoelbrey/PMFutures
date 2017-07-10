#!/bin/sh
#PBS -N Z3toP_2050RCP85
#PBS -l nodes=1:ppn=1
### Specify queue if necessary
#PBS -q batch
#PBS -j oe
#PBS -W group_list=fischer_group

cd /home/sbrey/projects/PMFutures/NCL

# Passed arguments: Variable, scenario, figure out how to pass these
# ncl  'scenario = "2000Firev1"' 'variable = "Z3"' 'fmod = "_fires_01."' hybridToPressure.ncl
# ncl  'scenario = "2000NoFirev1"' 'variable = "Z3"' 'fmod = "_nofires_01."' hybridToPressure.ncl
# ncl  'scenario = "2000FireHIv1"' 'variable = "Z3"' 'fmod = "_fires_hi_01."' hybridToPressure.ncl

ncl  'scenario = "2050RCP45Firev1"' 'variable = "Z3"' 'fmod = "_fires_01."' hybridToPressure.ncl

exit 0
