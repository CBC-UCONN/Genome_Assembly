#!/bin/bash
#SBATCH --job-name=polca
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=50G
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

OUTDIR=../../results/07_hornwort_hybrid/polish
mkdir -p $OUTDIR
cd $OUTDIR

GENOME=../masurca/CA.mr.99.17.15.0.02/primary.genome.scf.fasta
ILLUMINA=../../../data/illumina

#load software 

module load MaSuRCA/4.0.5
module load bwa/0.7.17

#run polca to polish

polca.sh \
	-a $GENOME \
	-r "$ILLUMINA/SRR10250248_1.fastq  $ILLUMINA/SRR10250248_2.fastq" \
	-t 32

