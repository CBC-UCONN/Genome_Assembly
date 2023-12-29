#!/bin/bash
#SBATCH --job-name=canu_10kb
#SBATCH -N 1               
#SBATCH -n 1               
#SBATCH -c 15           
#SBATCH --partition=general
#SBATCH --qos=general      
#SBATCH --mail-type=ALL    
#SBATCH --mem=35G          
#SBATCH --mail-user=
#SBATCH -o %x_%j.out                       
#SBATCH -e %x_%j.err 

module load gnuplot/5.2.2
module load canu/2.2

hostname
date

INDIR=../../results/02_Quality_Control/centrifuge
OUTDIR=../../results/03_Assembly/Canu_Assembly
mkdir -p $OUTDIR

canu -p canu -d $OUTDIR \
        genomeSize=600m \
        -minReadLength=10000 \
        -gridOptions="--partition=general --qos=general" canuIteration=1 \
        -nanopore $INDIR/filtered_Juglans_ailanthifolia.fasta  
