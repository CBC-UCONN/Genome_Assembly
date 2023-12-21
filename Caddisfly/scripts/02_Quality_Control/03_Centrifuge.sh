#!/bin/bash
#SBATCH --job-name=pacbio_centrifuge
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=100G
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load centrifuge/1.0.4-beta

#output directory
OUTDIR=../../results/02_Quality_Control/centrifuge
mkdir -p $OUTDIR

centrifuge -q \
  -p 12 \
  -x /isg/shared/databases/centrifuge/b+a+v+h/p_compressed+h+v \
  --report-file $OUTDIR/report.tsv \
  --quiet \
  --min-hitlen 50 \
        -k 1 \
  -U ../../data/SRR15654800.fastq.gz \
  --un $OUTDIR \
        >${OUTDIR}/classification.txt
