#!/bin/bash
#SBATCH --job-name=nanoplot
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo "HOSTNAME: `hostname`"
echo "Start Time: `date`"

module load NanoPlot/1.33.0

NanoPlot --fastq ../../data/nanopore/SRR9696346.fastq \
        --loglength \
        -o  ../../results/02_quality_control/nanopore \
        -t 10

date
module list
