#!/bin/bash
#SBATCH --job-name=QC_illumina_reads
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
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

INDIR=../../data/rnaseq

REPORTDIR=../../results/02_quality_control/rnaseq/fastp_reports
mkdir -p $REPORTDIR

TRIMDIR=../../results/02_quality_control/rnaseq/trimmed_sequences
mkdir -p $TRIMDIR

# run fastp to trip and generate reports on reads

fastp \
    --in1 $INDIR/SRR9662965_1.fastq.gz \
    --in2 $INDIR/SRR9662965_2.fastq.gz \
    --out1 $TRIMDIR/SRR9662965_trim_1.fastq.gz \
    --out2 $TRIMDIR/SRR9662965_trim_2.fastq.gz \
    --json $REPORTDIR/SRR9662965_fastp.json \
    --html $REPORTDIR/SRR9662965_fastp.html

module purge

########################################################
## Quality Control with fastqc 
#########################################################

module load fastqc/0.11.7

FASTQC=../../results/02_quality_control/rnaseq/fastqc_reports
mkdir -p $FASTQC

fastqc --outdir $FASTQC $TRIMDIR/SRR9662965_trim_{1..2}.fastq.gz
