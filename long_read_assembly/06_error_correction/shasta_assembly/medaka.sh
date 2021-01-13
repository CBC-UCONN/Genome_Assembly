#!/bin/bash
#SBATCH --job-name=medaka
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=180G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##              Medaka Assembly                         ##
##########################################################


module load medaka/0.11.4
module unload tabix/0.2.6

NPROC=16
BASECALLS=../../03_centrifuge/physcomitrellopsis_africana_rmv_contam.fasta
DRAFT=../../04_assembly/shasta/ShastaRun/Assembly.fasta
OUTDIR=$HOME/Genome_Assembly/long_read_assembly/06_error_correction/shasta_assembly

medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR}  -t 16


date
