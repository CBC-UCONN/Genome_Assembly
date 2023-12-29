#!/bin/bash
#SBATCH --job-name=centrifuge               
#SBATCH -N 1            
#SBATCH -n 1            
#SBATCH -c 12           
#SBATCH --partition=general      
#SBATCH --qos=general            
#SBATCH --mail-type=all          
#SBATCH --mem=100G               
#SBATCH --mail-user=
#SBATCH -o %x_%j.out                       
#SBATCH -e %x_%j.err

hostname
date

OUTDIR=../../results/02_Quality_Control/centrifuge
mkdir -p $OUTDIR

centrifuge -f \
        -x /core/labs/Wegrzyn/IngaGenome/Contam/longReads/f+b+a+v/abv \
        -p 12 \
        --report-file $OUTDIR/report.tsv \
        --quiet \
        --min-hitlen 50 \
        -U ../../data/Juglans_ailantifolia.fasta \
        --un $OUTDIR \
        >${OUTDIR}/classification.txt

date
