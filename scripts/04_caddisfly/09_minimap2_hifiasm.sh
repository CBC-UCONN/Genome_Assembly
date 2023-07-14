#!/bin/bash
#SBATCH --job-name=minimap_hifiasm
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
PACBIO=../../data/pacbio/SRR15654800.fastq.gz

OUTDIR=../../results/04_caddisfly/evaluation/hifiasm/alignment
    mkdir -p ${OUTDIR}

GENOME=../../results/04_caddisfly/hifiasm/hifiasm.fa

OUTROOT=SRR15654800


# run minimap
minimap2 -c --MD -ax map-hifi -t 15 ${GENOME} ${PACBIO} | \
samtools sort -@ 5 -T ${OUTDIR}/${OUTROOT}.temp -O BAM \
>${OUTDIR}/${OUTROOT}.bam

samtools index ${OUTDIR}/${OUTROOT}.bam

samtools stats ${OUTDIR}/${OUTROOT}.bam >${OUTDIR}/${OUTROOT}.stats
