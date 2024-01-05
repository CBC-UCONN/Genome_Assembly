#!/bin/bash
#SBATCH --job-name=medaka_quast
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


OUTDIR=../../results/07_Evaluation_Polish/quast
mkdir -p %OUTDIR

medaka=../../05_Error_Correction/medaka_flye_output/consensus.fasta



quast.py $medaka \
       --threads 10 \
       -o $OUTDIR


