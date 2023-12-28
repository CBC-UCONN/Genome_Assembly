#!/bin/bash
#SBATCH --job-name=walnut_symlink
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=your.email@uconn.edu
#SBATCH --mem=10G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

DATA=../../data/
mkdir -p $DATA


ln -s /archive/projects/PromethION/plant/JWalnut_assembly_files/Juglans_ailantifolia_pass.fastq.tar.gz $DATA/Juglans_ailantifolia_pass.fastq.tar.gz 



