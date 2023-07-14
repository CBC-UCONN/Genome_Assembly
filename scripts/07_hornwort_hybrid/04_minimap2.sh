#!/bin/bash
#SBATCH --job-name=minimap_ont
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
ONT=../../data/nanopore/SRR10190639_40.fastq

OUTDIR=../../results/07_hornwort_hybrid/masurca/alignment/minimap
    mkdir -p ${OUTDIR}

GENOME=../../results/07_hornwort_hybrid/masurca/CA.mr.99.17.15.0.02/primary.genome.scf.fasta

OUTROOT=SRR10190639_40_masurca


# run minimap
minimap2 -c --MD -ax map-ont -t 15 ${GENOME} ${ONT} | \
samtools sort -@ 5 -T ${OUTDIR}/${OUTROOT}.temp -O BAM \
>${OUTDIR}/${OUTROOT}.bam

samtools index ${OUTDIR}/${OUTROOT}.bam

samtools stats ${OUTDIR}/${OUTROOT}.bam >${OUTDIR}/${OUTROOT}.stats


