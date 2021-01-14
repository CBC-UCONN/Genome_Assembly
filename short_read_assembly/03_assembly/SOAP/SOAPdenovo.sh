#!/bin/bash
#SBATCH --job-name=SOAPdenovo
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


##########################################################
##		SOAP-denovo				##
##########################################################

module load SOAP-denovo/2.04

SOAPdenovo-127mer all \
        -s config_file \
        -K 131 \
        -p 8 \
        -R \
        -o graph_Sample_131 1>ass131.log 2>ass131.err

SOAPdenovo-127mer all \
        -s config_file \
        -K 135 \
        -p 8 \
        -R \
        -o graph_Sample_135 1>ass135.log 2>ass135.err
 
 
SOAPdenovo-127mer all \
        -s config_file \
        -K 141 \
        -p 8\
        -R \
        -o graph_Sample_141 1>ass141.log 2>ass141.err 

module unload SOAP-denovo/2.04

		
