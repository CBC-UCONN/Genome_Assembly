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

OUTDIR=../../../results/03_genome_size/caddisfly
    mkdir -p ${OUTDIR}
    cd ${OUTDIR}

INDIR=../../../data/pacbio

jellyfish count -t 30 -C -m 21 -s 100G -o $OUTDIR/21mer_out <(zcat $INDIR/SRR15654800.fastq.gz)
