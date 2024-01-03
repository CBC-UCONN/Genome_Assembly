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



OUTDIR=../../results/02_Quality_Control/kraken
mkdir -p $OUTDIR

INPUT=../../results/02_Quality_Control/trimmed_sequences

kraken2 -db /isg/shared/databases/kraken2/Standard \
	--paired $INPUT/SRR16841689_trim_1.fastq $INPUT/SRR16841689_trim_2.fastq \
	--use-names \
	--threads 16 \
	--output $OUTDIR/SRR16841689_kraken_general.out \
	--unclassified-out $OUTDIR/SRR16841689_unclassified#.fastq \
	--classified-out $OUTDIR/SRR16841689_classified#.fastq \
	--report $OUTDIR/SRR16841689_kraken_report.txt \
	--use-mpa-style 



date
