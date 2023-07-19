 # Genome Assembly   

This repository is a usable, publicly available tutorial. All steps have been provided for the UConn CBC Xanadu cluster here with appropriate headers for the Slurm scheduler that can be modified simply to run.  Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should use sbatch scriptname after modifying the script for each stage.  Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs.  If you are new to Linux, please use [this](https://bioinformatics.uconn.edu/unix-basics) handy guide for the operating system commands.  In this guide, you will be working with common bioinformatic file formats, such as [FASTA](https://en.wikipedia.org/wiki/FASTA_format), [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), [SAM/BAM](https://en.wikipedia.org/wiki/SAM_(file_format)), and [GFF3/GTF](https://en.wikipedia.org/wiki/General_feature_format). You can learn even more about each file format [here](https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/). If you do not have a Xanadu account and are an affiliate of UConn/UCHC, please apply for one **[here](https://bioinformatics.uconn.edu/contact-us/)**.   

**Contents**    

1.   [Data Download](#1-data-download)
2.   [Quality Control](#2-quality-control)
3.   [Genome Size Estimation](#3-genome-size-estimation)
4.   [Caddisfly Assembly](#4-caddisfly-assembly)
5.   [Hornwort Long Read Assembly](#5-hornwort-long-read-assembly)
6.   [Hornwort Short Read Assembly](#6-hornwort-short-read-assembly)
7.   [Hornwort Hybrid Assembly](#7-hornwort-hybrid-assembly)

## 1. Data Download

Follow the scripts in 01_data_download to obtain the publically available data. Here we are downloading pacbio HiFi, ONT, and illumina reads.


## 2. Quality Control

Follow the scripts in the 02_quality_control directory to characterize the reads with NanoPlot and FastQC. In here we also trim the illumina data and filter all the reads for contamination.


## 3. Genome Size Estimation

Follow the scripts in the 03_genome_size directory to generate counts of k-mers/ a histogram with Jellyfish to then take to GenomeScope. This will give an estimation for the genome size and heterozygosity of the sample.


## 4. Caddisfly Assembly

Follow the scripts in the 04_caddisfly directory to assemble a genome with PacBio HiFi data. We use both flye and hifiasm assemblers in this directory. We then evaluate with busco, quast, and minimap2.


## 5. Hornwort Long Read Assembly

In the 05_hornwort_ONT there are scripts to assemble a genome with ONT reads and Flye assembler. We then demonstrate polishing with medaka (long reads) and polca (short reads). We then evaluate with busco, quast, minimap2, and bwa.


## 6. Hornwort Short Read Assembly

In 06_hornwort_illumina we demonstrate a short read assembly with MaSuRCA. We then evaluate with busco, quast, and bwa.


## 7. Hornwort Hybrid Assembly

In 07_hornwort_hybrid we use both illumina and ONT reads to assemble a hybrid genome with MaSuRCA. We also demonstrate polishing with polca. We then evaluate with busco, quast, minimap2, and bwa.



