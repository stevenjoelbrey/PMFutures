#!/bin/sh
#PBS -N 2100RCP45_stag
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
# wind1000Lim = 8.   # m/s  
# wind500Lim  = 13.  # m/s
# precLim     = 0.01 # inches/day

# Passed arguments: hourlyVAR, scenario
python find_stagnation_days.py 2100RCP45 8. 13. 0.01

exit 0
