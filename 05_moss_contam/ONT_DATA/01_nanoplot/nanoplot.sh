#!/bin/bash
#SBATCH --job-name=nanoplot
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=50G
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo "HOSTNAME: `hostname`"
echo "Start Time: `date`"

module load


module load NanoPlot/1.21.0

NanoPlot --summary ../02_basecall_pass/sequencing_summary.txt \
        --loglength \
        -o summary-plots-log-transformed \
        -t 10

date
module list