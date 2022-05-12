#!/bin/bash
#SBATCH --job-name=01_bowtie2_align
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=firs.name@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

mkdir -p bowtie2_index

module load bowtie2/2.3.5.1
module load samtools/1.12

#bowtie2-build \
#	--threads 16 \
#	../02_masurca_assembly/CA.mr.41.17.20.0.02/final.genome.scf.fasta bowtie2_index/masurca_genome


#bowtie2 --threads 16 -x bowtie2_index/masurca_genome \
#	-1 ../01_ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L007_R1.fastq.gz,../01_ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L008_R1.fastq.gz \
#	-2 ../01_ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L007_R2.fastq.gz,../01_ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L008_R2.fastq.gz | \
#	samtools view -@ 16 -S -h -u - | \
#	samtools sort -@ 16 -T temp_bam > masurca_aligned.bam

samtools index -@ 16 masurca_aligned.bam masurca_aligned.bai




