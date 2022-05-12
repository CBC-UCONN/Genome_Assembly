#!/bin/bash
#SBATCH --job-name=02_kraken
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

module load kraken/2.0.8-beta
module load jellyfish/2.2.6

kraken2 -db /isg/shared/databases/kraken2/Standard \
	--gzip-compressed \
	--paired ../01_ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L007_R1.fastq.gz ../01_ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L007_R2.fastq.gz \
	--use-names \
	--threads 16 \
	--output S37_L007_kraken_general.out \
	--unclassified-out XUMNS_20180703_K00134_IL100105003_S37_L007_unclassified#.fastq \
	--classified-out XUMNS_20180703_K00134_IL100105003_S37_L007_classified#.fastq \
	--report S37_L007_kraken_report.txt \
	--use-mpa-style 

kraken2 -db /isg/shared/databases/kraken2/Standard \
	--gzip-compressed \
	--paired ../01_ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L008_R1.fastq.gz ../01_ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L008_R2.fastq.gz \
	--use-names \
        --threads 16 \
        --output S37_L008_kraken_general.out \
        --unclassified-out XUMNS_20180703_K00134_IL100105003_S37_L008_unclassified#.fastq \
        --classified-out XUMNS_20180703_K00134_IL100105003_S37_L008_classified#.fastq \
        --report S37_L008_kraken_report.txt \
        --use-mpa-style

date


