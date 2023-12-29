#!/bin/sh
#SBATCH --job-name=fastq_to_fasta
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH -c 1                         
#SBATCH --partition=general        
#SBATCH --qos=general                
#SBATCH --mail-type=END             
#SBATCH --mem=10G                   
#SBATCH --mail-user=stefan.wnuk@uconn.edu                
#SBATCH -o %x_%j.out              
#SBATCH -e %x_%j.err

module load seqkit/2.2.0

outdir=../../results/02_Quality_Control/centrifuge

seqtk seq -a Juglans_ailantifolia_filtered.fastq > filtered_Juglans_ailanthifolia.fasta
