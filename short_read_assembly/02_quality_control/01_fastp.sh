#!/bin/bash
#SBATCH --job-name=fastp_trimming
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
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

# set input/output directory variables
INDIR=../01_raw_reads/
REPORTDIR=fastp_reports
mkdir -p $REPORTDIR

# run fastp
fastp \
	--in1 $INDIR/SRR8776799_1.fastq \
	--in2 $INDIR/SRR8776799_2.fastq \
	--out1 trim_SRR8776799_1.fastq \
	--out2 trim_SRR8776799_2.fastq \
	--json $REPORTDIR/SRR8776799_fastp.json \
	--html $REPORTDIR/SRR8776799_fastp.html

