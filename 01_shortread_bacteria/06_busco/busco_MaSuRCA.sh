#!/bin/bash
#SBATCH --job-name=busco_MaSuRCA
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.name@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##		BUSCO 					##	
##########################################################
# MaSuRCA
module load busco/4.0.2


busco -i ../03_assembly/MaSuRCA/CA/final.genome.scf.fasta \
        -o MaSuRCA -l /isg/shared/databases/BUSCO/odb10/bacteria_odb10 -m genome


module load busco/4.0.2
