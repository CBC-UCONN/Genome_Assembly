#!/bin/bash
#SBATCH --job-name=busco_shasta
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
##              BUSCO                                   ##      
##########################################################

module load busco/4.0.2

export AUGUSTUS_CONFIG_PATH=../../augustus/config

busco -i ../../07_purge/shasta/curated.fasta \
        -o shasta_curated_busco -l /isg/shared/databases/BUSCO/odb10/viridiplantae_odb10 -m genome