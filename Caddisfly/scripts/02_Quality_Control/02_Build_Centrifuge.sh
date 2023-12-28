#!/bin/bash
#SBATCH --job-name=centrifuge_build
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=100G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date


OUTDIR=../../results/02_Quality_Control/Build_Centrifuge
mkdir -p $OUTDIR

cd $OUTDIR

module load centrifuge/1.0.4-beta
module load blast/2.13.0

centrifuge-download -o taxonomy taxonomy

centrifuge-download -o library -m -d "archaea,bacteria,viral" refseq > seqid2taxid.map

cat library/*/*.fna > input-sequences.fna

## build centrifuge index with 4 threads
centrifuge-build -p 4 --conversion-table seqid2taxid.map \
                 --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
                 input-sequences.fna abv
