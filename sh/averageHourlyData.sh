#!/bin/sh
#PBS -N 2100RCP45
#PBS -l nodes=1:ppn=1
### Specify queue if necessary
#PBS -q batch
#PBS -j oe
#PBS -W group_list=fischer_group

cd /home/sbrey/projects/PMFutures/Python

# Set the python Path so that it can use Anaconda distribution
export PATH="/home/sbrey/anaconda2/bin:$PATH"

# Passed arguments: hourlyVAR, scenario
python averageHourlyData.py PSL 2100RCP85

exit 0
