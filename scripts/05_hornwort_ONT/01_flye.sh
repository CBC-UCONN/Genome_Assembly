#!/bin/bash
#SBATCH --job-name=flye_nanopore
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=250G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##              Flye Assembly                           ##
##########################################################

module load flye/2.9.1

OUTDIR=../../results/05_hornwort_ONT/flye
mkdir -p $OUTDIR

flye --nano-raw ../../results/02_quality_control/nanopore/centrifuge/SRR10190639_40.fastq \
     --scaffold \
     --threads 32 \
     --out-dir $OUTDIR

