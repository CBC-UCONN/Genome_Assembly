#!/bin/bash
#SBATCH --job-name=quast_SPAdes
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=firs.name@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##		Quality Assesment: QUAST		##
##########################################################
# SPAdes
module load quast/5.0.2

quast.py ../03_assembly/SPAdes/scaffolds.fasta \
	--threads 8 \
	-o SPAdes


module unload quast/5.0.2

 

