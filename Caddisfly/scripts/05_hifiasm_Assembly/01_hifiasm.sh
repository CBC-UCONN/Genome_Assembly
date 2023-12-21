#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 30
#SBATCH --mem=500G
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

OUTDIR=../../results/05_hifiasm
mkdir -p $OUTDIR
cd $OUTDIR


PACBIO=../02_Quality_Control/centrifuge/pacbio_filtered.fastq

module load Hifiasm/0.18.5



hifiasm -o SRR15654800 -t 30 -l 2 $PACBIO
