#!/bin/bash
#SBATCH --job-name=nanopore_filter
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



##########################################
##   Extract Unclassified Sequences     ## 
##########################################

module load seqtk/1.2


OUTDIR=../../results/02_quality_control/nanopore/centrifuge
mkdir -p $OUTDIR
RAW=../../../../data/nanopore/SRR10190639_40.fastq


cd $OUTDIR

awk '($6/($7+1)) < 0.2' classification.txt | cut -f 1 >names.txt

seqtk subseq $RAW names.txt >SRR10190639_40.fastq


