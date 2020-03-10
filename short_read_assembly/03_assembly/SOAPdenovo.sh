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

SOAPdenovo-63mer all \
	-s config_file \
	-K 31 \
	-p 8 \
	-R \
	-o graph_Sample_31 1>ass31.log 2>ass31.err 

SOAPdenovo-63mer all \
	-s config_file \
	-K 35 \
	-p 8 \
	-R \
	-o graph_Sample_35 1>ass35.log 2>ass35.err 

SOAPdenovo-63mer all \
	-s config_file \
	-K 41 \
	-p 8\
	-R \
	-o graph_Sample_41 1>ass41.log 2>ass41.err 

module unload SOAP-denovo/2.04

		
