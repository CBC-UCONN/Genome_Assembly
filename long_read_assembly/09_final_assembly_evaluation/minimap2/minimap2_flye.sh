#!/bin/bash
#SBATCH --job-name=minimap2_flye
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

ref="../../08_purge/flye/curated.fasta"
read_file="../../03_centrifuge/physcomitrellopsis_africana_rmv_contam.fasta"

minimap2 -t 16 -ax map-ont ${ref} ${read_file} \
        | samtools sort -@ 16 -m 2G -o flye_aligned.bam -T flye_tmp.ali

module unload minimap2/2.15
module unload samtools/1.9

#Usage: bamtools stats [-in <filename> -in <filename> ... | -list <filelist>] [statsOptions]
module load bamtools/2.5.1
bamtools stats -in flye_aligned.bam