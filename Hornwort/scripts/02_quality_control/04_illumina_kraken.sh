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



OUTDIR=../../results/02_quality_control/illumina/kraken
mkdir -p $OUTDIR


kraken2 -db /isg/shared/databases/kraken2/Standard \
	--paired ../../results/02_quality_control/illumina/trimmed_sequences/SRR10250248_trim_1.fastq ../../results/02_quality_control/illumina/trimmed_sequences/SRR10250248_trim_2.fastq \
	--use-names \
	--threads 16 \
	--output $OUTDIR/SRR10250248_kraken_general.out \
	--unclassified-out $OUTDIR/SRR10250248_unclassified#.fastq \
	--classified-out $OUTDIR/SRR10250248_classified#.fastq \
	--report $OUTDIR/SRR10250248_kraken_report.txt \
	--use-mpa-style 



date
