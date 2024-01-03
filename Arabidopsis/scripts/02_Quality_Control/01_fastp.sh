#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`
date

#################################################################
# Trimming/QC of reads using fastp
#################################################################


module load fastp/0.23.2

INDIR=../../data

REPORTDIR=../../results/02_Quality_Control/fastp_reports
mkdir -p $REPORTDIR

TRIMDIR=../../results/02_Quality_Control/trimmed_sequences
mkdir -p $TRIMDIR

# run fastp to trim and generate reports on reads

fastp \
    --in1 $INDIR/SRR16841689_1.fastq \
    --in2 $INDIR/SRR16841689_2.fastq \
    --out1 $TRIMDIR/SRR16841689_trim_1.fastq \
    --out2 $TRIMDIR/SRR16841689_trim_2.fastq \
    --json $REPORTDIR/SRR16841689_fastp.json \
    --html $REPORTDIR/SRR16841689_fastp.html



