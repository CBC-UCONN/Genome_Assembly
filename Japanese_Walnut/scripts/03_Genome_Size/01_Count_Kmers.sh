#!/bin/bash
#SBATCH --job-name=kmer_filtered_genome_size
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=50G
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

module load gce/1.0.2

OUTDIR=../../results/03_Genome_Size_Estimation
mkdir -p $OUTDIR

cp read_files.lib $OUTDIR

cd $OUTDIR

# calculate k-mer frequency using kmerfreq
kmerfreq -k 17 -t 10 read_files.lib

# extract k-mer individual number
less read_files.lib.kmer.freq.stat | grep "#Kmer indivdual number"
less read_files.lib.kmer.freq.stat | perl -ne 'next if(/^#/ || /^\s/); print; ' | awk '{print $1"\t"$2}' > read_files.lib.kmer.freq.stat.2column 



