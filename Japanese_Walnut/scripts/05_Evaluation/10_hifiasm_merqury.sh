#!/bin/bash
#SBATCH --job-name=hifiasm_merqury
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

# Define input and output file names
filtered_kmer_db="../../results/05_Evaluation/meryl_db/kmer_db.filtered.meryl"
hifiasm="../../results/05_Evaluation/hifiasm/quast"
prefix="j_walnut_hifiasm"

# Execute merqury.sh
singularity exec /isg/shared/databases/nfx_singularity_cache/merqury.sif merqury.sh \
    $filtered_kmer_db \
    $hifiasm \
    $prefix
