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
mkdir -p SOAP_31_index SOAP_71_index SOAP_101_index

module load bowtie2/2.3.5.1

## SOAP_31
bowtie2-build \
	--threads 8 \
	../03_assembly/SOAP/graph_Sample_31.scafSeq SOAP_31_index/SOAP_31_index

bowtie2 -x SOAP_31_index/SOAP_31_index \
	-1 ../02_quality_control/trim_Sample_R1.fastq -2 ../02_quality_control/trim_Sample_R2.fastq \	
	-S SOAP_31.bowtie2.sam \
	--threads 8 2>SOAP_31.err

## SOAP_71
bowtie2-build \
        --threads 8 \
        ../03_assembly/SOAP/graph_Sample_71.scafSeq SOAP_71_index/SOAP_71_index

bowtie2 -x SOAP_71_index/SOAP_71_index \
        -1 ../02_quality_control/trim_Sample_R1.fastq -2 ../02_quality_control/trim_Sample_R2.fastq \
        -S SOAP_71.bowtie2.sam \
        --threads 8 2>SOAP_71.err

## SOAP_101 
bowtie2-build \
        --threads 8 \
	../03_assembly/SOAP/graph_Sample_101.scafSeq SOAP_101_index/SOAP_101_index

bowtie2 -x SOAP_101_index/SOAP_101_index \
	-1 ../02_quality_control/trim_Sample_R1.fastq -2 ../02_quality_control/trim_Sample_R2.fastq \
        -S SOAP_101.bowtie2.sam \
	--threads 8 2>SOAP_101.err



module unload bowtie2/2.3.5.1

 

