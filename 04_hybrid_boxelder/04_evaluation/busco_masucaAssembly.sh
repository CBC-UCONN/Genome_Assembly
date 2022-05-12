#!/bin/bash
#SBATCH --job-name=busco_masucaAssembly
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=20G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##              BUSCO                                   ##      
##########################################################

module load busco/5.0.0

export AUGUSTUS_CONFIG_PATH=../../augustus/config 

busco -i  ../02_masurca_assembly/CA.mr.41.17.20.0.02/final.genome.scf.fasta \
        -o busco_masucaAssembly -l /isg/shared/databases/BUSCO/odb10/lineages/viridiplantae_odb10 -m genome -c 16
