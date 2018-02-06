#!/bin/sh
#PBS -N ecoregion_1992
#PBS -l nodes=1:ppn=1
### Specify queue if necessary
#PBS -q batch
#PBS -j oe
#PBS -W group_list=fischer_group


cd /home/sbrey/projects/PMFutures

# Passed arguments: YYYY
Rscript --vanilla R/assign_ecoregion_to_FPAFOD.R 1992

exit 0