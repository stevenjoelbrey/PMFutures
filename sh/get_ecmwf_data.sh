#!/bin/sh
#PBS -N get_era_1996
#PBS -l nodes=1:ppn=1
### Specify queue if necessary
#PBS -q batch
#PBS -j oe
#PBS -W group_list=fischer_group


cd /home/sbrey/projects/PMFutures/Python

# Set the python Path so that it can use Anaconda distribution
export PATH="/home/sbrey/anaconda2/bin:$PATH"


# Passed arguments: hourlyVAR, scenario
python get_ERA_Interim_data.py z 1996 1996 pl

exit 0
