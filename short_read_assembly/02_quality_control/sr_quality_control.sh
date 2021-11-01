#!/bin/bash
#SBATCH --job-name=sr_quality_control
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo "HOSTNAME: `hostname`"
echo "Start Time: `date`"


##########################################################
## quality control                                      ##
##########################################################

module load Trimmomatic/0.39

java -jar $Trimmomatic PE -threads 4 \
        ../01_raw_reads/Sample_R1.fastq \
        ../01_raw_reads/Sample_R2.fastq \
        trim_Sample_R1.fastq trim_Sample_R1_singles.fastq \
        trim_Sample_R2.fastq trim_Sample_R2_singles.fastq \
        ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 \
        SLIDINGWINDOW:4:25 MINLEN:45


module unload Trimmomatic/0.39
echo "End of trimming: `date`"


##########################################################
## FASTQC Raw Reads                                     ##
##########################################################
module load fastqc/0.11.7
mkdir -p RAWfastqc_OUT

fastqc -o ./RAWfastqc_OUT ../01_raw_reads/Sample_R1.fastq ../01_raw_reads/Sample_R2.fastq


##########################################################
## FASTQC Trimmed Reads                                 ##
##########################################################
mkdir -p TRIMfastqc_OUT

fastqc -o ./TRIMfastqc_OUT ./trim_Sample_R1.fastq ./trim_Sample_R2.fastq


echo "End of FASTQC: `date`"
