#!/bin/bash
#SBATCH -J ragtag
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=100G
#SBATCH --mail-type=END
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date
module load RagTag/2.1.0

reference=/core/projects/EBP/conservation/can_butternut/juglans_reference/GCA_022457165.1_NFU_Jman_1.0_genomic.fna
purged_canu=../../results/08_Purge_Haplotigs/Canu/curated.fasta
OUTDIR=../../results/10_Scaffolding
mkdir -p $OUTDIR

# scaffold a query assembly
ragtag.py scaffold \
	$reference \
	$purged_canu \
	-o $OUTDIR \
	-f 1000
