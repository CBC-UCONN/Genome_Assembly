#!/bin/bash
#SBATCH --job-name=busco_SOAP
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=firs.name@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##		BUSCO					##	
##########################################################
# SOAP

module load busco/4.0.2 

#export AUGUSTUS_CONFIG_PATH=/isg/shared/apps/augustus/3.2.3/config 

busco -i ../03_assembly/SOAP/graph_Sample_31.scafSeq \
        -o SOAP_31 -l /isg/shared/databases/BUSCO/odb10/bacteria_odb10 -m genome


busco -i ../03_assembly/SOAP/graph_Sample_71.scafSeq \
        -o SOAP_71 -l /isg/shared/databases/BUSCO/odb10/bacteria_odb10 -m genome


busco -i ../03_assembly/SOAP/graph_Sample_101.scafSeq \
        -o SOAP_101 -l /isg/shared/databases/BUSCO/odb10/bacteria_odb10 -m genome

module unload busco/4.0.2

date

