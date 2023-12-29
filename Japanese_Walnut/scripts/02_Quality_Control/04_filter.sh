#!/bin/bash
#SBATCH --job-name=filter
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=20G
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

outdir=../../results/02_Quality_Control/centrifuge

cd $outdir

module load seqkit/2.2.0

#taking the classified reads from the STDOUT file and adding those to a new file 
grep -vw "unclassified" classification.txt > contaminated_reads.txt


#extracting just the IDs from the text file
awk NF=1 contaminated_reads.txt > contaminated_read_ids.txt


#keeping only unique IDs and removing any duplicates
sort -u contaminated_read_ids.txt > no_dup_centrifuge_contaminated_read_ids.txt


#adds only the sequences that DO NOT match the headers in the file we provided
seqkit grep -v -f no_dup_centrifuge_contaminated_read_ids.txt ../../../data/Juglans_ailantifolia.fastq.gz  > Juglans_ailantifolia_filtered.fastq  

