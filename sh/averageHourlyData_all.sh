#!/bin/bash



for NSR in 2000Base 2050RCP45 2050RCP85 2100RCP45 2100RCP85
do
 cat > makeAverageHourlyData.sh <<EOF                                                                           
#!/bin/bash                                                                                                        
#PBS -N MURPH$NSR
### Number of nodes and processors per node.                                                                        
### This job will run on 1 nodes with 1 cores per node                                                             
#PBS -l nodes=1:ppn=1                                                                                               
### Specify queue if necessary
##PBS -q batch                                                                                                      
### Specified resources (walltime, memory, etc) if necessary                                                        
### Join standard output and standard error in one file                                                             
#PBS -j oe
#PBS -W group_list=fischer_group


echo started at `date`
export PATH="/home/sbrey/anaconda2/bin:$PATH"
cd /home/sbrey/projects/PMFutures/Python
ulimit -s unlimited                                                                                                 

python averageHourlyData.py PSL $NSR 2000 

echo ended at `date`
 
EOF

 qsub makeAverageHourlyData.sh

done
