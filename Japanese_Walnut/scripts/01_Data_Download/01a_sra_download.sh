#!/bin/bash
#SBATCH --job-name=caddisfly_download
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 6
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=your.email@uconn.edu
#SBATCH --mem=40G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date
