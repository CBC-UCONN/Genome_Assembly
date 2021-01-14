#!/bin/bash
#SBATCH --job-name=quast_SOAP
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=firs.name@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##		Quality Assesment: QUAST		##
##########################################################

module load quast/5.0.2

quast.py ../03_assembly/SOAP/graph_Sample_31.scafSeq \
        --threads 8 \
        -o SOAP_31

quast.py ../03_assembly/SOAP/graph_Sample_71.scafSeq \
        --threads 8 \
        -o SOAP_71

quast.py ../03_assembly/SOAP/graph_Sample_101.scafSeq \
        --threads 8 \
        -o SOAP_101



module unload quast/5.0.2

 

