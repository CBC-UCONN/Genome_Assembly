#!/bin/bash
#SBATCH --job-name=Kraken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date


cat ../01_Trimmomatic/573_paired_1.fq ../01_Trimmomatic/574_paired_1.fq ../01_Trimmomatic/575_paired_1.fq ../01_Trimmomatic/576_paired_1.fq ../01_Trimmomatic/577_paired_1.fq ../01_Trimmomatic/578_paired_1.fq ../01_Trimmomatic/579_paired_1.fq ../01_Trimmomatic/580_paired_1.fq >> cat_R1.fq

cat ../01_Trimmomatic/573_paired_2.fq ../01_Trimmomatic/574_paired_2.fq ../01_Trimmomatic/575_paired_2.fq ../01_Trimmomatic/576_paired_2.fq ../01_Trimmomatic/577_paired_2.fq ../01_Trimmomatic/578_paired_2.fq ../01_Trimmomatic/579_paired_2.fq ../01_Trimmomatic/580_paired_2.fq >> cat_R2.fq


module load kraken/2.0.8-beta
module load jellyfish/2.2.6

kraken2 -db /isg/shared/databases/kraken2/Standard \
        --paired cat_R1.fq cat_R2.fq \
        --use-names \
        --threads 16 \
        --output kraken_general.out \
        --unclassified-out unclassified#.fastq \
        --classified-out classified#.fastq      \
        --report kraken_report.txt \
        --use-mpa-style

date
