#!/bin/bash
#SBATCH --job-name=busco
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=20G
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

module load busco/4.0.2

export AUGUSTUS_CONFIG_PATH=../../augustus/config 

busco -i  ../05_purge/curated.fasta \
	-o busco_purge -l /isg/shared/databases/BUSCO/odb10/lineages/viridiplantae_odb10 -m genome -c 8


