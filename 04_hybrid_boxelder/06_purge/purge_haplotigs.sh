#!/bin/bash
#SBATCH --job-name=02purge_haplotigs
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=40G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

##########################################
##      minimap alignment               ##
##########################################
#module load minimap2/2.15
#module load samtools/1.9

ref="../05_polish/pilon_out/masurca_polished.fasta"
read_file="../01_Pacbio_DATA/acne_pb.fasta"

align_file="masurca_aligned"

#https://bitbucket.org/mroachawri/purge_haplotigs/src/master/
#minimap2 -t 16 -ax map-pb ${ref} ${read_file} \
#	| samtools sort -@ 16 -m 1G -o ${align_file}.bam -T tmp.ali

#samtools faidx ${ref} 
#date
#samtools index ${align_file}.bam ${align_file}.bai
#date
#module unload minimap2/2.15
#module unload samtools/1.9


##########################################
##      purge haplotigs                 ## 
##########################################
module load R/4.1.2 
module load purge_haplotigs/1.0
## STEP-1
## purge_haplotigs  readhist  -b aligned.bam  -g genome.fasta  [ -t threads ]
#purge_haplotigs readhist -b ${align_file}.bam -g ${ref} -t 16

## STEP-2
## purge_haplotigs  contigcov  -i aligned.bam.genecov  -l 30  -m 80  -h 145  [-o coverage_stats.csv -j 80  -s 80 ]
#purge_haplotigs  contigcov -i ${align_file}.bam.gencov -l 20 -m 75 -h 190  

## STEP-3
purge_haplotigs purge -b ${align_file}.bam -g ${ref} -c coverage_stats.csv -d -a 60


