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
mkdir -p SOAP_31_index SOAP_35_index SOAP_41_index

module load bowtie2/2.3.5.1

## SOAP_31
bowtie2-build \
	--threads 8 \
	../03_assembly/SOAP/graph_Sample_31.scafSeq SOAP_31_index/SOAP_31_index

bowtie2 -x SOAP_31_index/SOAP_31_index \
	-U ../01_raw_reads/Sample_R1.fastq,../01_raw_reads/Sample_R2.fastq \
	-S SOAP_31.bowtie2.sam \
	--threads 8 2>SOAP_31.err

## SOAP_35
bowtie2-build \
        --threads 8 \
        ../03_assembly/SOAP/graph_Sample_35.scafSeq SOAP_35_index/SOAP_35_index

bowtie2 -x SOAP_35_index/SOAP_35_index \
        -U ../01_raw_reads/Sample_R1.fastq,../01_raw_reads/Sample_R2.fastq \
        -S SOAP_35.bowtie2.sam \
        --threads 8 2>SOAP_35.err

## SOAP_41 
bowtie2-build \
        --threads 8 \
	../03_assembly/SOAP/graph_Sample_41.scafSeq SOAP_41_index/SOAP_41_index

bowtie2 -x SOAP_41_index/SOAP_41_index \
	-U ../01_raw_reads/Sample_R1.fastq,../01_raw_reads/Sample_R2.fastq \
        -S SOAP_41.bowtie2.sam \
	--threads 8 2>SOAP_41.err



module unload bowtie2/2.3.5.1

 

