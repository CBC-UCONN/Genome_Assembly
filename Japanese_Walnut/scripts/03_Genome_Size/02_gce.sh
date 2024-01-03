#!/bin/bash
#SBATCH --job-name=gce
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=50G
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load gce/1.0.2


OUTDIR=../../results/03_Genome_Size_Estimation
mkdir -p $OUTDIR

cd $OUTDIR


# run GCE in homozygous mode
gce -g 122579822915 -f read_files.lib.kmer.freq.stat.2column >gce.table 2>gce.log

