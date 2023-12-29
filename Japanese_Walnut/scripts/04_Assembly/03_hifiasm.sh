#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mem=250G
#SBATCH -o %x_%j.out                       
#SBATCH -e %x_%j.err 


module load Hifiasm/0.18.5

reads=../../results/02_Quality_Control/centrifuge/Juglans_ailantifolia_filtered.fastq
OUTDIR=../../results/03_Assembly/hifiasm

hifiasm -o $OUTDIR/japanese_walnut -t 32 $reads
