#!/bin/bash
#SBATCH --job-name=extract_fasta
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date


##convert gfa output to fa


WORKDIR=../../results/05_hifiasm

cd $WORKDIR

awk '/^S/{print ">"$2;print $3}' SRR15654800.bp.p_ctg.gfa > hifiasm.fa
