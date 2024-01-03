#!/bin/bash
#SBATCH --job-name=merqury
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=128G
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=END
#SBATCH --mail-user=
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

# Load Singularity
module load singularity

OUTDIR=../../results/07_Evaluation_Polish/merqury 
mkdir -p $OUTDIR

cd $OUTDIR

# Define input and output file names
filtered_kmer_db="../../05_Evaluation/meryl_db/kmer_db.filtered.meryl"
flye="../../05_Error_Correction/medaka_flye_output/consensus.fasta"
prefix="flye_medaka_merqury"

# Execute merqury.sh
singularity exec /isg/shared/databases/nfx_singularity_cache/merqury.sif merqury.sh \
    $filtered_kmer_db \
    $flye \
    $prefix
