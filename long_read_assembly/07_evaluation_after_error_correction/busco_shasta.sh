#!/bin/bash
#SBATCH --job-name=busco_shasta
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --partition=general
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

module load busco/4.0.2

export AUGUSTUS_CONFIG_PATH=../../augustus/config

busco -i ../../06_error_correction/shasta_assembly/consensus.fasta \
        -o busco_shasta -l /isg/shared/databases/BUSCO/odb10/viridiplantae_odb10 -m genome