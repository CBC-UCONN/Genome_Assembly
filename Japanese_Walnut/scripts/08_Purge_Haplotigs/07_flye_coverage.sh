#!/bin/bash
#BATCH --job-name=flye_contigcov
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

ref=../../04_Assembly/Flye/assembly.fasta 
WORKDIR=../../results/08_Purge_Haplotigs/Flye 

cd $WORKDIR


##########################################
##      purge haplotigs                 ##
##########################################
module load purge_haplotigs/1.1.2
module load R/4.2.2

## Step 2: use the cutoff values from histogram 
purge_haplotigs  contigcov -i flye_aligned.bam.gencov -l 40 -m 160 -h 288


