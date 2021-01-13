#!/bin/bash
#SBATCH --job-name=bowtie2_SOAP
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=firs.name@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##		Quality Assesment: bowtie2		##
##########################################################
# SOAP
mkdir -p SOAP_131_index SOAP_135_index SOAP_141_index

module load bowtie2/2.3.5.1

## SOAP_131
bowtie2-build \
        --threads 8 \
        ../03_assembly/SOAP/graph_Sample_131.scafSeq SOAP_131_index/SOAP_131_index

bowtie2 -x SOAP_131_index/SOAP_131_index \
        -1 ../01_raw_reads/Sample_R1.fastq -2 ../01_raw_reads/Sample_R2.fastq \
        -S SOAP_131.bowtie2.sam \
        --threads 8 2>SOAP_131.err

## SOAP_135
bowtie2-build \
        --threads 8 \
        ../03_assembly/SOAP/graph_Sample_135.scafSeq SOAP_135_index/SOAP_135_index

bowtie2 -x SOAP_135_index/SOAP_135_index \
        -1 ../01_raw_reads/Sample_R1.fastq -2 ../01_raw_reads/Sample_R2.fastq \
        -S SOAP_135.bowtie2.sam \
        --threads 8 2>SOAP_135.err

## SOAP_141 
bowtie2-build \
        --threads 8 \
        ../03_assembly/SOAP/graph_Sample_141.scafSeq SOAP_141_index/SOAP_141_index

bowtie2 -x SOAP_141_index/SOAP_141_index \
        -1 ../01_raw_reads/Sample_R1.fastq -2 ../01_raw_reads/Sample_R2.fastq \
        -S SOAP_141.bowtie2.sam \
        --threads 8 2>SOAP_141.err



module unload bowtie2/2.3.5.1

 

