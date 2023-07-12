#!/bin/bash
#SBATCH --job-name=sra_download
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


# data
# Illumina (Illumina HiSeq 1500) Genome sequence of Anthoceros agrestis Oxford strain: 400-bp insert PCR-Free Discovar protocol
    # bioproject: PRJNA574453
        # biosample: SAMN12858332
        # SRA runs:
            # SRR10250248
            
# Nanopore (MiniION) sequencing of Anthoceros agrestis Oxford
    # bioproject: PRJNA574453
        # biosample: SAMN12858332
        # SRA runs:
            # SRR10190639
            # SRR10190640
            


# load software
module load sratoolkit/3.0.5

# output directory, create if it doesn't exist

OUTDIR=../../data/illumina
mkdir -p ${OUTDIR}

cd ${OUTDIR}

fasterq-dump SRR10250248


# nanopore outdir

NANODIR=../nanopore
mkdir -p ${NANODIR}

cd ${NANODIR}

fasterq-dump SRR10190639
fasterq-dump SRR10190640
