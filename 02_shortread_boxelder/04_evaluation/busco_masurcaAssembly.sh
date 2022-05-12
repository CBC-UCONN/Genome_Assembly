#!/bin/bash
#SBATCH --job-name=busco_masurcaAssembly
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=30G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=firs.name@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##              BUSCO                                   ##      
##########################################################

module load busco/5.0.0 

export AUGUSTUS_CONFIG_PATH=../../augustus/config 

busco -i  ../03_assembly/02_masurca/CA/final.genome.scf.fasta \
	-c 16 \
	-o busco_masurcaAssembly -l /isg/shared/databases/BUSCO/odb10/lineages/viridiplantae_odb10 -m genome 


