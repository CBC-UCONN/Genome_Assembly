#!/bin/bash
#BATCH --job-name=purge_haplotigs
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

ref="../../06_error_correction/flye_assembly/consensus.fasta"
read_file="../../03_centrifuge/physcomitrellopsis_africana_rmv_contam.fasta"

minimap2 -t 16 -ax map-ont ${ref} ${read_file} \
       | samtools view -hF 256 - \
       | samtools sort -@ 16 -m 1G -o flye_aligned.bam -T flye_tmp.ali 

module unload minimap2/2.15
module unload samtools/1.9

##########################################
##      purge haplotigs                 ## 
##########################################
module load purge_haplotigs/1.0
## STEP-1
## purge_haplotigs  readhist  -b aligned.bam  -g genome.fasta  [ -t threads ]
purge_haplotigs readhist -b flye_aligned.bam -g ${ref} -t 16

## STEP-2
## purge_haplotigs  contigcov  -i aligned.bam.genecov  -l 30  -m 80  -h 145  [-o coverage_stats.csv -j 80  -s 80 ]
purge_haplotigs  contigcov -i flye_aligned.bam.gencov -l 3 -m 55 -h 195  

## STEP-3
purge_haplotigs purge -b flye_aligned.bam -g ${ref} -c coverage_stats.csv -d -a 60