#!/bin/sh
#SBATCH --job-name=fastq_to_fasta
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH -c 1                         
#SBATCH --partition=general        
#SBATCH --qos=general                
#SBATCH --mail-type=END             
#SBATCH --mem=10G                   
#SBATCH --mail-user=               
#SBATCH -o %x_%j.out              
#SBATCH -e %x_%j.err

hostname
date

data=../../data
RAW=../../data/Juglans_ailantifolia_pass.fastq.tar.gz

tar -xvzf $RAW

mv Juglans_ailantifolia_pass.fastq.gz $data/Juglans_ailantifolia.fastq.gz
FASTQ=../../data/Juglans_ailantifolia.fastq.gz

module load seqtk/1.3

seqtk seq -a $FASTQ > $data/Juglans_ailantifolia.fasta
