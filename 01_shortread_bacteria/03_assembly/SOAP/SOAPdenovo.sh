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
        -K 31 \
        -p 8 \
        -R \
        -o graph_Sample_31 1>ass31.log 2>ass31.err 

SOAPdenovo-127mer all \
        -s config_file \
        -K 71 \
        -p 8 \
        -R \
        -o graph_Sample_71 1>ass71.log 2>ass71.err 

SOAPdenovo-127mer all \
        -s config_file \
        -K 101 \
        -p 8\
        -R \
        -o graph_Sample_101 1>ass101.log 2>ass101.err 

module unload SOAP-denovo/2.04

		
