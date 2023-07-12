#!/bin/bash
#SBATCH --job-name=medaka
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=80G
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# software
module load medaka/1.7.3

# input/output
BASECALLS=../../data/nanopore/SRR10190639_40.fastq
DRAFT=../../results/05_hornwort_ONT/flye/assembly.fasta
OUTDIR=../../results/05_hornwort_ONT/medaka

NPROC=16
MODEL=r941_prom_sup_g507 

# run guppy
medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t ${NPROC} -m ${MODEL}
