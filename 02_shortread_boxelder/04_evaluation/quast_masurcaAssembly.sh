#!/bin/bash
#SBATCH --job-name=quast_masurcaAssembly
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=20G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=neranjan.perera@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##              QUAST                                   ##      
##########################################################

module load quast/5.0.2

quast.py ../03_assembly/02_masurca/CA/final.genome.scf.fasta \
        --threads 8 \
        -o quast_masurcaAssembly


