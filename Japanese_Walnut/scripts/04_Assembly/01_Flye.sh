#!/bin/bash
#SBATCH --job-name=flye
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=300G
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=END
#SBATCH --mail-user=anthony.he@uconn.edu
#SBATCH -o %x_%j.out                       
#SBATCH -e %x_%j.err 

module load flye/2.9.1

INDIR=../../results/02_Quality_Control/centrifuge
OUTDIR=../../results/04_Assembly/Flye
flye --nano-hq $INDIR/filtered_Juglans_ailanthifolia.fasta \
  --no-alt-contigs \
  --threads 32 \
  --out-dir $OUTDIR \
  --scaffold \
  --asm-coverage 60 \
  --genome-size 572M
