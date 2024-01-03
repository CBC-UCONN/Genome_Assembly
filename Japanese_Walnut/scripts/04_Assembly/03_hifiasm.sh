#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=150G
#SBATCH -o %x_%j.out                       
#SBATCH -e %x_%j.err 


module load Hifiasm/0.18.5

reads=../../results/02_Quality_Control/centrifuge/Juglans_ailantifolia_filtered.fastq
OUTDIR=../../results/03_Assembly/hifiasm
mkdir -p $OUTDIR

hifiasm -o $OUTDIR/japanese_walnut -t 32 $reads
