#!/bin/bash
#BATCH --job-name=flye_purge
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

WORKDIR=../../results/08_Purge_Haplotigs/Flye 

cd $WORKDIR

## Step 3: final step that purges seqs based on your cutoffs 
purge_haplotigs purge -b flye_aligned.bam -g ${ref} -c coverage_stats.csv -d
