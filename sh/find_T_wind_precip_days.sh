#!/bin/sh
#PBS -N 2100RCP85_TWP
#PBS -l nodes=1:ppn=1
### Specify queue if necessary
#PBS -q batch
#PBS -j oe
#PBS -W group_list=fischer_group


cd /home/sbrey/projects/PMFutures/Python

# Set the python Path so that it can use Anaconda distribution
export PATH="/home/sbrey/anaconda2/bin:$PATH"


# Passed arguments: hourlyVAR, scenario
python find_T_wind_precip_days.py useValue 2100RCP85 298 9 0.01

exit 0
