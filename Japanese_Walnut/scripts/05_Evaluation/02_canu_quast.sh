#!/bin/bash
#SBATCH --job-name=canu_quast
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=15G
#SBATCH --mail-user=
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

module load quast/5.2.0


CANU=../../results/05_Evaluation/Canu/quast
mkdir -p $CANU


quast.py ../../results/04_Assembly/Canu/canu.contigs.fasta \
       --threads 10 \
       -o $CANU


