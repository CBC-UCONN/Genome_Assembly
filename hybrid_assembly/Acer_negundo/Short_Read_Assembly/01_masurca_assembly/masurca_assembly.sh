#!/bin/bash
#SBATCH --job-name=masurca_assembly
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=150G
#SBATCH --mail-user=neranjan.perera@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date


module load singularity/3.1.1
module load MaSuRCA/3.5.0
module load perl/5.30.1

masurca config.txt

./assemble.sh


