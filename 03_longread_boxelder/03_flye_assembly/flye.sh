#!/bin/bash
#SBATCH --job-name=flye_pacbio
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=400G
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=neranjan.perera@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##              Flye Assembly                           ##
##########################################################

module load flye/2.4.2

flye --pacbio-raw ../01_Pacbio_DATA/acne_pb.fasta \
        --genome-size 300m \
        --threads 32 \
        --out-dir . 

