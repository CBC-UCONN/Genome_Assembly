#!/bin/bash
#SBATCH --job-name=quast_flye
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

GENOME=../../results/05_hornwort_ONT/flye/assembly.fasta
OUTDIR=../../results/05_hornwort_ONT/evaluation/quast
    mkdir -p ${OUTDIR}

quast.py ${GENOME} \
        --threads 16 \
        -o ${OUTDIR}
