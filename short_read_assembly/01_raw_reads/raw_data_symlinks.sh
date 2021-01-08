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

fpath="/UCHC/PublicShare/CBC_Tutorials/Genome_Assembly/short_read_assembly/01_raw_reads/"

for f in ${fpath}*; do
        echo $f
        echo `basename ${f}`
        ln -s ${f} `basename ${f}`
done
