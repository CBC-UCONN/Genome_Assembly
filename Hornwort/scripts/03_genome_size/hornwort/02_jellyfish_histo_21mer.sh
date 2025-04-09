#!/bin/bash
#SBATCH --job-name=21_mer_histo_hornwort
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 30
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# software
module load jellyfish/2.2.6

# files/directories
INDIR=../../../results/03_genome_size/hornwort

# run jellyfish
jellyfish histo -o ${INDIR}/21mer.histo ${INDIR}/21mer_out
