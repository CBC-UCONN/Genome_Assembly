#!/bin/bash
#SBATCH --job-name=kraken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

module load kraken/2.1.2
module load jellyfish/2.3.0

kraken2 -db /isg/shared/databases/kraken2/Standard \
	--gzip-compressed \
	--paired ../../results/02_quality_control/rnaseq/trimmed_sequences/SRR9662965_trim_1.fastq.gz ../../results/02_quality_control/rnaseq/trimmed_sequences/SRR9662965_trim_2.fastq.gz \
	--use-names \
	--threads 16 \
	--output SRR9662965_kraken_general.out \
	--unclassified-out SRR9662965_unclassified#.fastq \
	--classified-out SRR9662965_classified#.fastq \
	--report SRR9662965_kraken_report.txt \
	--use-mpa-style 



date
