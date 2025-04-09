#!/bin/bash
#SBATCH --job-name=busco_polca
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##              BUSCO                                   ##      
##########################################################

module load busco/5.4.5

OUTDIR="../../results/07_hornwort_hybrid/evaluation/polca/busco"
    mkdir -p ${OUTDIR}
GENOME=../../results/07_hornwort_hybrid/polca/primary.genome.scf.fasta.PolcaCorrected.fa
DATABASE="/isg/shared/databases/BUSCO/odb10/lineages/viridiplantae_odb10"

busco \
    -i ${GENOME} \
    -o ${OUTDIR} \
    -l ${DATABASE} \
    -m genome \
    -c 8 \
    -f

