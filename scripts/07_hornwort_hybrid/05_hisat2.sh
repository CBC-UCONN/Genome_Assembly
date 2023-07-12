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
OUTDIR=../../results/07_hornwort_hybrid/alignments/index
mkdir -p $OUTDIR

GENOME=../../results/07_hornwort_hybrid/masurca/CA.mr.99.17.15.0.02/primary.genome.scf.fasta
ILLUMINA=../../results/02_quality_control/illumina/trimmed_sequences


hisat2-build -p 16 $GENOME $OUTDIR/Aagr

##Next Align

INDEX=../../results/07_hornwort_hybrid/alignments/index/Aagr


ALIGNMENTS=../../results/07_hornwort_hybrid/alignments
hisat2 \
    -x $INDEX \
    -1 $ILLUMINA/SRR10250248_trim_1.fastq \
    -2 $ILLUMINA/SRR10250248_trim_2.fastq | \
samtools view -@ 1 -S -h -u - | \
samtools sort -@ 1 -T SRR10250248 - >$ALIGNMENTS/SRR10250248.bam

# index bam files
samtools index $ALIGNMENTS/SRR10250248.bam

#stats
samtools stats $ALIGNMENTS/SRR10250248.bam >$ALIGNMENTS/SRR10250248.stats
