#!/bin/bash
#SBATCH --job-name=bowtie2_MaSuRCA
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=firs.name@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##		Quality Assesment: bowtie2		##
##########################################################
# MaSuRCA
mkdir -p MaSuRCA_index

module load bowtie2/2.3.5.1

bowtie2-build \
	--threads 8 \
	../03_assembly/MaSuRCA/CA/final.genome.scf.fasta MaSuRCA_index/MaSuRCA_index

bowtie2 -x MaSuRCA_index/MaSuRCA_index \
        -1 ../01_raw_reads/Sample_R1.fastq -2 ../01_raw_reads/Sample_R2.fastq \
        -S MaSuRCA.bowtie2.sam \
        --threads 8 2>MaSuRCA.err


module unload bowtie2/2.3.5.1

 

