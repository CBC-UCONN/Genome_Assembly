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

INDIR=../../data/illumina

REPORTDIR=../../results/02_quality_control/illumina/fastp_reports
mkdir -p $REPORTDIR

TRIMDIR=../../results/02_quality_control/illumina/trimmed_sequences
mkdir -p $TRIMDIR

# run fastp to trim and generate reports on reads

fastp \
    --in1 $INDIR/SRR10250248_1.fastq \
    --in2 $INDIR/SRR10250248_2.fastq \
    --out1 $TRIMDIR/SRR10250248_trim_1.fastq \
    --out2 $TRIMDIR/SRR10250248_trim_2.fastq \
    --json $REPORTDIR/SRR10250248_fastp.json \
    --html $REPORTDIR/SRR10250248_fastp.html

module purge

########################################################
## Quality Control with fastqc 
#########################################################

module load fastqc/0.11.7

FASTQC=../../results/02_quality_control/illumina/fastqc_reports
mkdir -p $FASTQC

fastqc --outdir $FASTQC $TRIMDIR/SRR10250248_trim_{1..2}.fastq
