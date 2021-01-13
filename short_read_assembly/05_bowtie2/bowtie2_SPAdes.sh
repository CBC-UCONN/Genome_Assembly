#!/bin/bash
#SBATCH --job-name=bowtie2_SPAdes
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
##		Quality Assesment: bowtie2		##
##########################################################
# SPAdes
mkdir -p SPAdes_index

module load bowtie2/2.3.5.1

bowtie2-build \
	--threads 8 \
	../03_assembly/SPAdes/scaffolds.fasta SPAdes_index/SPAdes_index

bowtie2 -x SPAdes_index/SPAdes_index \
        -1 ../01_raw_reads/Sample_R1.fastq -2 ../01_raw_reads/Sample_R2.fastq \
        -S SPAdes.bowtie2.sam \
        --threads 8 2>SPAdes.err

module unload bowtie2/2.3.5.1

 

