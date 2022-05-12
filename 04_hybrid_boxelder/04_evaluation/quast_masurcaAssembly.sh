#!/bin/bash
#SBATCH --job-name=masurcaAssembly
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=20G
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

quast.py ../02_masurca_assembly/CA.mr.41.17.20.0.02/final.genome.scf.fasta \
        --threads 16 \
        -o quast_masurcaAssembly
