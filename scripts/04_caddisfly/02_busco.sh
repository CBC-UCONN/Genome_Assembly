#!/bin/bash
#SBATCH --job-name=busco_flye
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

OUTDIR="../../results/04_caddisfly/evaluation/flye/busco"
    mkdir -p ${OUTDIR}
GENOME="../../results/04_caddisfly/flye/assembly.fasta"
DATABASE="/isg/shared/databases/BUSCO/odb10/lineages/insecta_odb10"

busco \
    -i ${GENOME} \
    -o ${OUTDIR} \
    -l ${DATABASE} \
    -m genome \
    -c 8 \
    -f


