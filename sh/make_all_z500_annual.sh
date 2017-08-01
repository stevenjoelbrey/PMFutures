#!/bin/sh
#PBS -N regrid_fire_2016
#PBS -l nodes=1:ppn=1
### Specify queue if necessary
#PBS -q batch
#PBS -j oe
#PBS -W group_list=fischer_group

cd /home/sbrey/projects/PMFutures/Python

# Set the python Path so that it can use Anaconda distribution
export PATH="/home/sbrey/anaconda2/bin:$PATH"

# Passed arguments: hourlyVAR, scenario
# "year", "species"
python make_all_z500_annual.py 

exit 0
