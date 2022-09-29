#!/bin/bash
#SBATCH --job-name=SPAdes
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=firs.name@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##		SPAdes					##
##########################################################

module load SPAdes/3.13.0

spades.py \
	-1 ../../02_quality_control/trim_SRR8776799_1.fastq \
	-2 ../../02_quality_control/trim_SRR8776799_2.fastq \
	--careful \
	--threads 8 \
	--memory 30 \
	-o . 


module unload SPAdes/3.13.0

