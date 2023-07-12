#!/bin/bash
#SBATCH --job-name=masurca_assembly
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=150G
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

OUTDIR=../../results/06_hornwort_illumina/masurca
mkdir -p $OUTDIR
cp config.txt $OUTDIR
cd $OUTDIR

module load MaSuRCA/4.0.9

masurca config.txt

./assemble.sh
