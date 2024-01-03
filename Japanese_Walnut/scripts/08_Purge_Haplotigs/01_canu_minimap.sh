#!/bin/bash
#BATCH --job-name=canu_minimap
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=200G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

module load minimap2/2.15
module load samtools/1.9

ref=../../results/04_Assembly/Canu/canu.contigs.fasta
read_file=../../results/02_Quality_Control/centrifuge/filtered_Juglans_ailantifolia.fasta
OUTDIR=../../results/08_Purge_Haplotigs/Canu
mkdir -p $OUTDIR

minimap2 -t 16 -ax map-ont ${ref} ${read_file} \
       | samtools view -hF 256 - \
       | samtools sort -@ 16 -m 1G -o $OUTDIR/canu_aligned.bam -T canu_tmp.ali



