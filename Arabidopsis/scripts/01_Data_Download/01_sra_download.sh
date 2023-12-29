#!/bin/bash
#SBATCH --job-name=Arabidopsis_download
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 6
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=your.email@uconn.edu
#SBATCH --mem=40G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# run this script to download the data from the SRA


# Arabidopsis thaliana: young true leaves
    # bioproject: PRJNA777107


# load software
module load sratoolkit/3.0.2

# output directory, create if it doesn't exist

OUTDIR=../../data
mkdir -p $OUTDIR

cd ${OUTDIR}

#PacBio HIFI 
	#SRA: SRR18777598

fastq-dump --gzip SRR18777598
mv SRR18777598.fastq.gz Athal_pacbio.fastq.gz


#ONT AB 3130xL Genetic Analyzer
	#SRA: SRR16832054

fasterq-dump SRR16832054


# Illumina Genome Analyzer IIx
	#SRA: SRR16841689

fasterq-dump SRR16841689

date
