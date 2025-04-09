#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

#######index genome#######

module load bwa/0.7.17
module load samtools/1.16.1


INDEXDIR=../../results/05_hornwort_ONT/evaluation/flye_raw/alignment/bwa/index
mkdir -p $INDEXDIR
GENOME=../../results/05_hornwort_ONT/flye/assembly.fasta

bwa index \
	-p $INDEXDIR/Aagr \
	$GENOME


##################################
# align sequences to the reference
##################################

# align sequences to the reference genome

# input, output files and directories
FQ1=../../results/02_quality_control/illumina/trimmed_sequences/SRR10250248_trim_1.fastq
FQ2=../../results/02_quality_control/illumina/trimmed_sequences/SRR10250248_trim_2.fastq

OUTDIR=../../results/05_hornwort_ONT/evaluation/flye_raw/alignment/bwa
mkdir -p $OUTDIR


# get sample ID
SAM=SRR10250248


# align sequences--
# run bwa mem to align, then pipe it to samtools

bwa mem -t 4 $INDEXDIR/Aagr $FQ1 $FQ2 | \
	samtools sort -@ 5 -T ${OUTDIR}/${SAM}.temp -O BAM \
	>${OUTDIR}/${SAM}.bam

samtools index ${OUTDIR}/${SAM}.bam

samtools stats ${OUTDIR}/${SAM}.bam >${OUTDIR}/${SAM}.stats
