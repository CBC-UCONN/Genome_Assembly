#!/bin/bash
#SBATCH --job-name=canu_busco
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 15
#SBATCH --mem=20G
#SBATCH --partition=xeon #pay attention!!! we are using xeon, not general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=keertana.chagari@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err



module load busco/5.4.5

busco -i /core/projects/EBP/conservation/japanese_walnut/Student_Japanese_walnut_dir/assembly/07_purge_haplotigs/Canu_Purge/curated.fasta \
        -o canu_busco_output -l embryophyta_odb10 -m genome -c 15
