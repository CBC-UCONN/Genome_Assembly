#!/bin/bash
#SBATCH --job-name=flye_pacbio
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

OUTDIR=../../results/04_Flye
mkdir -p $OUTDIR

flye --pacbio-hifi ../../results/02_Quality_Control/centrifuge/pacbio_filtered.fastq \
     --scaffold \
     --genome-size 350m \
     --threads 32 \
     --out-dir $OUTDIR
