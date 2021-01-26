#!/bin/bash
#SBATCH --job-name=quast
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
##              QUAST                                   ##      
##########################################################

module load quast/5.0.2

quast.py ../01_masurca_assembly/CA/final.genome.scf.fasta \
        --threads 16 \
        -o masurca_assembly