#!/bin/sh
#PBS -N make_daily_vo
#PBS -l nodes=1:ppn=1
### Specify queue if necessary
#PBS -q batch
#PBS -j oe
#PBS -W group_list=fischer_group


cd /home/sbrey/projects/PMFutures/Python

# Set the python Path so that it can use Anaconda distribution
export PATH="/home/sbrey/anaconda2/bin:$PATH"


# Passed arguments: hourlyVAR, scenario
python average6HourlyData.py vo 1997 2016

exit 0
