#!/bin/bash
#SBATCH --job-name=01_qc
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo "HOSTNAME: `hostname`"
echo "Start Time: `date`"

##########################################################
## FASTQC Raw Reads                                     ##
##########################################################
module load fastqc/0.11.7
mkdir -p RAWfastqc_OUT

fastqc -t 2 -o ./RAWfastqc_OUT ../01_ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L007_R1.fastq.gz ../01_ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L007_R2.fastq.gz 

fastqc -t 2 -o ./RAWfastqc_OUT ../01_ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L008_R1.fastq.gz ../01_ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L008_R2.fastq.gz

module unload fastqc/0.11.7

##########################################################
## MultiQC
##########################################################
module load MultiQC/1.10.1
multiqc --title RawMultiQC --filename RawMultiQC --outdir RawMultiQC ./RAWfastqc_OUT

date



