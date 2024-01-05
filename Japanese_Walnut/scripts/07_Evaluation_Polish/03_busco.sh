#!/bin/bash
#SBATCH --job-name=medaka_busco
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 15
#SBATCH --mem=20G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

module load quast/5.2.0

OUTDIR=../../results/07_Evaluation_Polish/busco
mkdir -p $OUTDIR

GENOME=../../05_Error_Correction/medaka_flye_output/consensus.fasta \

module load busco/5.4.5

busco -i $GENOME \
      -o $OUTDIR \
      -l embryophyta_odb10 \
      -m genome \
      -c 15
