#!/bin/bash
#SBATCH --job-name=testflye_pacbio
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

OUTDIR=../../results/04_caddisfly/test
mkdir -p $OUTDIR

flye --pacbio-hifi ../../data/pacbio/SRR15654800.fastq.gz \
     --scaffold \
     --genome-size 350m \
     --threads 32 \
     --out-dir $OUTDIR

