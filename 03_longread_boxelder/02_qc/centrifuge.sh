#!/bin/bash
#SBATCH --job-name=centrifuge
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
#SBATCH --mail-type=END
#SBATCH --mem=100G
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load centrifuge/1.0.4-beta

centrifuge -f \
	-p 12 \
	-x /isg/shared/databases/centrifuge/b+a+v+h/p_compressed+h+v \
	--report-file report.tsv \
	--quiet \
	--min-hitlen 50 \
	-U ../01_Pacbio_DATA/acne_pb.fasta \
	--un . \
	--un-conc .



date
 
