#!/bin/bash
#SBATCH --job-name=TESTmasurca_assembly
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=100G
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

OUTDIR=../../results/07_hornwort_hybrid/masurca_test
mkdir -p $OUTDIR
cp config.txt $OUTDIR
cd $OUTDIR


module load MaSuRCA/4.1.0


masurca config.txt

./assemble.sh
