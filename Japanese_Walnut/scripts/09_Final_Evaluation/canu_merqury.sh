#!/bin/bash
#SBATCH --job-name=canu_merqury
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=128G
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=END
#SBATCH --mail-user=keertana.chagari@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

# Load Singularity
module load singularity

# Define input and output file names
output_filtered_kmer_db="/core/projects/EBP/conservation/japanese_walnut/Student_Japanese_walnut_dir/assembly/04_inicial_assembly_evaluation/meryl_db/kmer_db.filtered.meryl/"
output_flye="/core/projects/EBP/conservation/japanese_walnut/Student_Japanese_walnut_dir/assembly/07_purge_haplotigs/Canu_Purge/curated.fasta"
output_prefix="canu_purge"

# Execute merqury.sh
singularity exec /isg/shared/databases/nfx_singularity_cache/merqury.sif merqury.sh \
    $output_filtered_kmer_db \
    $output_flye \
    $output_prefix
