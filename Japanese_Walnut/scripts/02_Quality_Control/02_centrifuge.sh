#!/bin/bash
#SBATCH --job-name=centrifuge               
#SBATCH -N 1            
#SBATCH -n 1            
#SBATCH -c 10           
#SBATCH --partition=general      
#SBATCH --qos=general            
#SBATCH --mail-type=all          
#SBATCH --mem=100G               
#SBATCH --mail-user=
#SBATCH -o %x_%j.out                       
#SBATCH -e %x_%j.err       


hostname
date

module load centrifuge/1.0.4-beta

# Run the centrifuge program with the following parameters:
# -f: specify the input file format as FASTA
# -x: specify the centrifuge index database to use
# --report-file: specify the filename of the output report
# --quiet: suppress non-error messages
# --min-hitlen: specify the minimum length of a hit to be considered
# -U: specify the input file (FASTA)

OUTDIR=../../results/02_Quality_Control/centrifuge

centrifuge -f \
        -x /core/labs/Wegrzyn/IngaGenome/Contam/longReads/f+b+a+v/abv \
        --report-file Centrifugereport.tsv \
        --quiet \
        --min-hitlen 50 \
        -U ../../data/Juglans_ailantifolia.fasta



