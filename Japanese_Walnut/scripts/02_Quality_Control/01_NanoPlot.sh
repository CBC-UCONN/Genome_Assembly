#!/bin/bash
#SBATCH --job-name=NanoPlot
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=80G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo "HOSTNAME: `hostname`"
echo "Start Time: `date`"


module load NanoPlot/1.33.0

DATA=../../data

OUTDIR=../../results/02_Quality_Control/NanoPlot

NanoPlot --fastq $DATA/Juglans_ailantifolia.fastq.gz  --loglength --verbose -o $OUTDIR/summary-fastqpass -t 10 -p summary-fastqpass 


# Generate a summary plot of the raw reads with log-transformed data
NanoPlot --fasta $DATA/Juglans_ailantifolia.fasta -o $OUTDIR/summary-txt -t 10


echo "End Time: `date`"


