#!/bin/bash
#SBATCH --job-name=medaka
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=180G
#SBATCH --partition=xeon
#SBATCH --qos=general                          
#SBATCH --mail-type=ALL                      
#SBATCH --mail-user=
#SBATCH -o %x_%A.out                         
#SBATCH -e %x_%A.err                         

module load medaka/1.7.1

OUTDIR=../../results/06_Error_Correction
mkdir -p $OUTDIR

medaka_consensus -t 32 -i ../../results/02_Quality_Control/centrifuge/filtered_Juglans_ailantifolia.fasta -d ../../results/04_Assembly/Flye/assembly.fasta -o $OUTDIR/medaka

