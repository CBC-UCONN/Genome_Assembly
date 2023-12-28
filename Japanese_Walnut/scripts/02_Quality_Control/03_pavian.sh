#!/bin/bash
#SBATCH --job-name=Japanese_walnut_pavian_report
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=10G
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

module load centrifuge/1.0.4-beta

index=/core/labs/Wegrzyn/IngaGenome/Contam/longReads/f+b+a+v

centrifuge-kreport -x $index/abv ../../results/02_Quality_Control/centrifuge/Centrifugereport.tsv > pavian_report

