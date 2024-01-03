#!/bin/bash
#SBATCH --job-name=meryl
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=128G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

# Load Singularity
module load singularity

# Define input and output file names
input="../../data/Juglans_ailantifolia_pass_filtered.fastq.gz"
output_kmer_db="../../results/05_Evaluation/meryl_db/kmer_db.meryl"
output_filtered_kmer_db="../../results/05_Evaluation/meryl_db/kmer_db.filtered.meryl"

# Execute meryl count
singularity exec /isg/shared/databases/nfx_singularity_cache/merqury.sif meryl count \
    threads=32 \
    k=19 \
    $input \
    output $output_kmer_db

# Execute meryl greater-than
singularity exec /isg/shared/databases/nfx_singularity_cache/merqury.sif meryl greater-than 1 \
    threads=32 \
    k=19 \
    output $output_filtered_kmer_db $output_kmer_db
