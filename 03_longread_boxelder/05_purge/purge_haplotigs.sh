#!/bin/bash
#SBATCH --job-name=04purge_haplotigs
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
module load minimap2/2.15
module load samtools/1.9

ref="../03_flye_assembly/scaffolds.fasta"
read_file="../01_Pacbio_DATA/acne_pb.fasta"

#https://bitbucket.org/mroachawri/purge_haplotigs/src/master/
minimap2 -t 16 -ax map-pb ${ref} ${read_file} \
	| samtools sort -@ 16 -m 1G -o flye_aligned.bam -T flye_tmp.ali

samtools faidx ${ref} 

samtools index flye_aligned.bam flye_aligned.bai

module unload minimap2/2.15
module unload samtools/1.9


##########################################
##      purge haplotigs                 ## 
##########################################
module load R/4.1.2 
module load purge_haplotigs/1.0
## STEP-1
## purge_haplotigs  readhist  -b aligned.bam  -g genome.fasta  [ -t threads ]
purge_haplotigs readhist -b flye_aligned.bam -g ${ref} -t 16

## STEP-2
## purge_haplotigs  contigcov  -i aligned.bam.genecov  -l 30  -m 80  -h 145  [-o coverage_stats.csv -j 80  -s 80 ]
purge_haplotigs  contigcov -i flye_aligned.bam.gencov -l 15 -m 70 -h 190  

## STEP-3
purge_haplotigs purge -b flye_aligned.bam -g ${ref} -c coverage_stats.csv -d -a 60


