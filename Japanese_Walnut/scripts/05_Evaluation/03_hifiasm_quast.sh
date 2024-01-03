#!/bin/bash
#SBATCH --job-name=hifiasm_quast
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

HIFI=../../results/05_Evaluation/hifiasm/quast
mkdir -p $HIFI

quast.py ../../results/04_Assembly/hifiasm/japanese_walnut.asm.bp.p_ctg.fa \
       --threads 10 \
       -o $HIFI


