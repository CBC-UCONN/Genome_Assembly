#!/bin/bash
#SBATCH --job-name=caddisfly_download
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


# PacBio HiFi genome of Atopsyche davidsoni: larvae whole body
    # bioproject: PRJNA741212
        # biosample: SAMN20492399 
        # SRA runs:
            # SRR15654800

# load software
module load sratoolkit/3.0.2

# output directory, create if it doesn't exist


# pacbio outdir

DATA=../../data
mkdir -p ${DATA}

cd ${DATA}

fastq-dump --gzip SRR15654800


date

