#!/bin/bash
#SBATCH --job-name=fastqer_dump
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=2G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

#################################################################
# Download fastq files from SRA 
#################################################################

module load sratoolkit/2.11.3

# Staphylococcus aureus
    # BioProject: PRJNA528186
    # BioSample: SAMN11175590
    # Illumina data SRA run ID: SRR8776799 

# download SRA run
fasterq-dump SRR8776799

