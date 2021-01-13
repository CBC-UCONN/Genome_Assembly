#!/bin/bash
#SBATCH --job-name=centrifuge
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=END
#SBATCH --mem=200G
#SBATCH --mail-user=neranjan.perera@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

module load centrifuge/1.0.4-beta

centrifuge -f \
        -x /isg/shared/databases/centrifuge/b+a+v+h/p_compressed+h+v \
        --report-file report.tsv \
        --quiet \
        --min-hitlen 50 \
        -U ../02_basecall_pass/5074_test_LSK109_30JAN19-reads-pass.fasta

date
