#!/bin/bash
#SBATCH --job-name=02_pilon
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=300G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=firs.name@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date


module load pilon/1.24

export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch/${USER}

mkdir -p pilon_out 

java -Xmx200g -XshowSettings -jar $PILON \
	--genome ../02_masurca_assembly/CA.mr.41.17.20.0.02/final.genome.scf.fasta \
	--frags masurca_aligned.bam \
	--output masurca_polished \
	--outdir pilon_out



