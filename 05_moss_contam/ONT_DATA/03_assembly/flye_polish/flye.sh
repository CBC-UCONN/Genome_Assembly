#!/bin/bash
#SBATCH --job-name=flye
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 64
#SBATCH --mem=500G
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##              Flye Assembly                           ##
##########################################################

module load flye/2.4.2

flye --nano-raw ../../03_centrifuge/physcomitrellopsis_africana_rmv_contam.fasta \
        --genome-size 500m \
        --threads 32 \
        --iterations 3 \
        --out-dir .
 

#--iterations 3 : 3 round of polishing