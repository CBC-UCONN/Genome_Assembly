#!/bin/sh
#SBATCH --job-name=fastq_to_fasta
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH -c 1                         
#SBATCH --partition=general        
#SBATCH --qos=general                
#SBATCH --mail-type=END             
#SBATCH --mem=50G                   
#SBATCH --mail-user=stefan.wnuk@uconn.edu                
#SBATCH -o %x_%j.out              
#SBATCH -e %x_%j.err

hostname
date

module load seqtk/1.3

outdir=../../results/02_Quality_Control/centrifuge

cd $outdir

seqtk seq -a Juglans_ailantifolia_filtered.fastq > filtered_Juglans_ailanthifolia.fasta
