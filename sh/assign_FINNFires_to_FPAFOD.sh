#!/bin/sh
#PBS -N match_FINN_fires_2007
#PBS -l nodes=1:ppn=12
### Specify queue if necessary
#PBS -q batch
#PBS -j oe
#PBS -W group_list=fischer_group


cd /home/sbrey/projects/PMFutures

Rscript --vanilla R/assign_FINNFires_to_FPAFOD.R  2007

exit 0