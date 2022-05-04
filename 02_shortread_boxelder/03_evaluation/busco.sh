#!/bin/bash
#SBATCH --job-name=busco
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.name@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##              BUSCO                                   ##      
##########################################################

module load busco/4.0.2

export AUGUSTUS_CONFIG_PATH=/UCHC/PublicShare/CBC_Tutorials/Genome_Assembly/long_read_assembly/augustus/config 

busco -i  ../01_masurca_assembly/CA/final.genome.scf.fasta \
        -o busco -l /isg/shared/databases/BUSCO/odb10/viridiplantae_odb10 -m genome 