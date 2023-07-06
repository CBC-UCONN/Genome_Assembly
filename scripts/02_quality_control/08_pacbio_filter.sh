#!/bin/bash
#SBATCH --job-name=pacbio_filter
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


OUTDIR=../../results/02_quality_control/pacbio/centrifuge
mkdir -p $OUTDIR
RAW=../../../../data/pacbio/SRR15654800.fastq.gz

cd $OUTDIR

awk '$6 < 5000' classification.txt | cut -f 1 >names.txt

seqtk subseq $RAW names.txt >pacbio_filtered.fastq

chmod +x pacbio_filtered.fastq
