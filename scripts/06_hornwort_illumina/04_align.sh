#!/bin/bash
#SBATCH --job-name=hisat2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=80G
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load hisat2/2.2.1
module load samtools/1.16.1

# input/output directories
OUTDIR=../../results/06_hornwort_illumina/alignments/index
mkdir -p $OUTDIR

GENOME=../../results/06_hornwort_illumina/masurca/CA/primary.genome.scf.fasta
ILLUMINA=../../results/02_quality_control/illumina/trimmed_sequences


hisat2-build -p 16 $GENOME $OUTDIR/Aagr

##Next Align

INDEX=../../results/06_hornwort_illumina/alignments/index/Aagr

SAMPLE=SRR10250248

ALIGNMENTS=../../results/06_hornwort_illumina/alignments
hisat2 \
    -x $INDEX \
    -1 $ILLUMINA/${SAMPLE}_trim_1.fastq \
    -2 $ILLUMINA/${SAMPLE}_trim_2.fastq | \
samtools view -@ 1 -S -h -u - | \
samtools sort -@ 1 -T $SAMPLE - >$ALIGNMENTS/$SAMPLE.bam

# index bam files
samtools index $ALIGNMENTS/$SAMPLE.bam

#stats
samtools stats $ALIGNMENTS/${SAMPLE}.bam >$ALIGNMENTS/${SAMPLE}.stats
