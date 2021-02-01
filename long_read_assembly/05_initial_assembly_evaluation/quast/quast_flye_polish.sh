#!/bin/bash
#SBATCH --job-name=quast_flye_polish
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
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
##              QUAST                                   ##      
##########################################################

module load quast/5.0.2

quast.py ../../04_assembly/flye_polish/assembly.fasta \
        --threads 16 \
        -o quast_flye_polish