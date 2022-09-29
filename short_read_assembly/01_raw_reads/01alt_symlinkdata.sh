#!/bin/bash
#SBATCH --job-name=raw_data_symlinks
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=128M
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# If you're working on Xanadu, you should 
# run this script to symlink the data from 
# our raw data directory rather than downloading
# it from the SRA

# Staphylococcus aureus
    # BioProject: PRJNA528186
    # BioSample: SAMN11175590
    # Illumina data SRA run ID: SRR8776799 

fpath="/core/cbc/tutorials/rawdata/Genome_Assembly/staph_aureus/"

for f in ${fpath}*; do
        echo $f
        echo `basename ${f}`
        ln -s ${f} `basename ${f}`
done
