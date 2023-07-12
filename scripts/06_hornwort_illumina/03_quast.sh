#!/bin/bash
#SBATCH --job-name=quast
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=20G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##              QUAST                                   ##      
##########################################################

module load quast/5.2.0

GENOME=../../results/06_hornwort_illumina/masurca/CA/primary.genome.scf.fasta
OUTDIR=../../results/06_hornwort_illumina/evaluation/quast
    mkdir -p ${OUTDIR}

quast.py ${GENOME} \
        --threads 16 \
        -o ${OUTDIR}


