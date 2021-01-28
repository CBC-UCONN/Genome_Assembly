#!/bin/bash
#SBATCH --job-name=masurca_hybrid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=ALL
#SBATCH --mem=500G
#SBATCH --mail-user=neranjan.perera@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

module load singularity/3.1.1
module load MaSuRCA/3.5.0

module load perl/5.30.1

masurca config.txt

./assemble.sh
