#!/bin/bash
#SBATCH --job-name=flye_quast
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=15G
#SBATCH --mail-user=sefan.wnuk@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

module load quast/5.2.0

quast.py --threads 10 -o flye_quast_output /core/projects/EBP/conservation/japanese_walnut/Student_Japanese_walnut_dir/assembly/07_purge_haplotigs/Flye_Purge/curated.fasta

