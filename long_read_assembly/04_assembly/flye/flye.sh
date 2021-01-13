#!/bin/bash
#SBATCH --job-name=flye
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=600G
#SBATCH --partition=himem2
#SBATCH --qos=himem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=firs.name@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##              Flye Assembly                           ##
##########################################################

module load flye/2.4.2

flye --nano-raw ../../03_centrifuge/physcomitrellopsis_africana_rmv_contam.fasta \
        --genome-size 1g \
        --threads 32 \
        --out-dir /UCHC/PublicShare/CBC_Tutorials/Genome_Assembly/long_read_assembly/03_assembly/flye_t