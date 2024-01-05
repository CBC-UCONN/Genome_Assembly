#!/bin/bash
#BATCH --job-name=flye_readhist
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=200G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

ref=../../04_Assembly/Flye/assembly.fasta 
WORKDIR=../../results/08_Purge_Haplotigs/Flye 

cd $WORKDIR

##########################################
##      purge haplotigs                 ##
##########################################
module load purge_haplotigs/1.1.2
module load R/4.2.2

## STEP-1: you just created a histogram; look at it and find high, med, and low cutoffs for next step :)
purge_haplotigs readhist -b flye_aligned.bam -g ${ref} -t 16 -d 300


