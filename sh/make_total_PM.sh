#!/bin/sh
#PBS -N 2000Firev1
#PBS -l nodes=1:ppn=1
### Specify queue if necessary
#PBS -q batch
#PBS -j oe
#PBS -W group_list=fischer_group

cd /home/sbrey/projects/PMFutures/Python

# Set the python Path so that it can use Anaconda distribution
export PATH="/home/sbrey/anaconda2/bin:$PATH"

# Argument order 
# scenario    = '2000Base'

# Passed arguments: hourlyVAR, scenario
python make_total_PM.py 2000Firev1

exit 0
