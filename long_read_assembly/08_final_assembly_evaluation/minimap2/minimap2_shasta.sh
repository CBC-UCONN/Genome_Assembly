#!/bin/bash
#SBATCH --job-name=minimap2_shasta
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=40G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=neranjan.perera@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

##########################################
##      minimap alignment               ##
##########################################
module load minimap2/2.15
module load samtools/1.9

ref="../../07_purge/shasta/curated.fasta"
read_file="../../02_basecall_pass/5074_test_LSK109_30JAN19-reads-pass.fasta"

minimap2 -t 16 -ax map-ont ${ref} ${read_file} \
        | samtools sort -@ 16 -m 2G -o shasta_aligned.bam -T shasta_tmp.ali

module unload minimap2/2.15
module unload samtools/1.9



#Usage: bamtools stats [-in <filename> -in <filename> ... | -list <filelist>] [statsOptions]
module load bamtools/2.5.1
bamtools stats -in shasta_aligned.bam