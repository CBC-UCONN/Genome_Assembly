#!/bin/bash
#SBATCH --job-name=shasta
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=500G
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.name@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##              shasta Assembly                         ##
##########################################################

module load shasta/0.5.1

shasta --input ../../02_basecall_pass/5074_test_LSK109_30JAN19-reads-pass.fasta \
        --Reads.minReadLength 500 \
        --memoryMode anonymous \
        --threads 32

date
