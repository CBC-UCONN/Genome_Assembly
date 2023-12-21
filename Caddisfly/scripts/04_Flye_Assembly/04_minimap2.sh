#!/bin/bash
#SBATCH --job-name=minimap
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 15
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=80G
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load software

module load minimap2/2.24
module load samtools/1.16.1

# input/output files, directories
PACBIO=../../results/02_Quality_Control/centrifuge/pacbio_filtered.fastq

OUTDIR=../../results/04_Flye/alignment
    mkdir -p ${OUTDIR}

GENOME=../../results/04_Flye/assembly.fasta

OUTROOT=SRR15654800


# run minimap
minimap2 -c --MD -ax map-hifi -t 15 ${GENOME} ${PACBIO} | \
samtools sort -@ 5 -T ${OUTDIR}/${OUTROOT}.temp -O BAM \
>${OUTDIR}/${OUTROOT}.bam

samtools index ${OUTDIR}/${OUTROOT}.bam

samtools stats ${OUTDIR}/${OUTROOT}.bam >${OUTDIR}/${OUTROOT}.stats
