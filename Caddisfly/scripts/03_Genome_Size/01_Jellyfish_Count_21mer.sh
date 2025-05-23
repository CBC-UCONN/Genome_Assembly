#!/bin/bash
#SBATCH --job-name=jellyfish_pacbio
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 30
#SBATCH --mem=400G
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load software
module load jellyfish/2.2.6

# files/directories

OUTDIR=../../results/03_Genome_Size
    mkdir -p ${OUTDIR}

INDIR=../../results/02_Quality_Control/centrifuge

jellyfish count -t 30 -C -m 21 -s 100G -o $OUTDIR/21mer_out $INDIR/pacbio_filtered.fastq

# -C, --canonical, count both strands
# -m, --mer-size, size of the k-mer
