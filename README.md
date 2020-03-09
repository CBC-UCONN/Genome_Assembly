# Genome_Assembly   

This repository is a usable, publicly available tutorial. All steps have been provided for the UConn CBC Xanadu cluster here with appropriate headers for the Slurm scheduler that can be modified simply to run.  Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should use sbatch scriptname after modifying the script for each stage.  Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs.  If you are new to Linux, please use [this](https://bioinformatics.uconn.edu/unix-basics) handy guide for the operating system commands.  In this guide, you will be working with common bioinformatic file formats, such as [FASTA](https://en.wikipedia.org/wiki/FASTA_format), [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), [SAM/BAM](https://en.wikipedia.org/wiki/SAM_(file_format)), and [GFF3/GTF](https://en.wikipedia.org/wiki/General_feature_format). You can learn even more about each file format [here](https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/). If you do not have a Xanadu account and are an affiliate of UConn/UCHC, please apply for one **[here](https://bioinformatics.uconn.edu/contact-us/)**.   

**Contents**    

1.   [Overview](#1-overview)  
2.   [Short Read Genome Assembly](#2-short-read-genome-assembly)   


## 1. Overview   
This tutorial will teach you how to use open source quality control, genome assembly, and assembly assessment tools to complete a high quality de novo assembly which is commonly utilized when you dont have a reference genome. Moving through the tutorial, you will take pair end short read data from a bacterial species and perform assemblies via various commonly used genome assmeblers. With these assemblies completed we will then need to assess the quality of the data produced. Once finished with the short read data we will move on to performing a long read assembly with long read nanopore data using basecalling and commonly used long read assemblers. Finally, we will then move on to Hybrid PacBio data.     

![](images/overview-1.png)

#### Structure of the tutorial   

The tutorial is organized into 3 parts:   
*   Short Read Assembly   
*   Long Read Assembly  
*   Hybrid Assembly   

Once you git clone the repository, you can see three main folders:  
```
Genome_Assembly/
├── short_read_assembly/
├── long_read_assembly/
└── hybrid_assembly/ 
```

#### SLURM scripts

 In each folder it will contain the scripts to run each job in the HPC cluster. The tutorial will be using SLURM schedular to submit jobs to Xanadu cluster. In each script we will be using it will contain a header section which will allocate the resources for the SLURM schedular. The header section will contain:   
```
#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
```
Before beginning, we need to understand a few aspects of the Xanadu server. When first logging into Xanadu from your local terminal, you will be connected to the submit node. The submit node is the interface with which users on Xanadu may submit their processes to the desired compute nodes, which will run the process. Never, under any circumstance, run processes directly in the submit node. Your process will be killed and all of your work lost! This tutorial will not teach you shell script configuration to submit your tasks on Xanadu. Therefore, before moving on, read and master the topics covered in the [Xanadu tutorial](https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/xanadu/).  


 ## 2. Short Read Genome Assembly   
 When using short reads for genome assembly, the reads first need to be quality controled, and then will be assembled. In this section we will be showing how to use three assemblers: [SOAPdenovo](https://github.com/aquaskyline/SOAPdenovo2), [SPAdes](https://github.com/ablab/spades) and [MaSuRCA](https://github.com/alekseyzimin/masurca) to assemble short reads and then how to evaluate the assembled reads.   

 ![](images/short_read_assembly-1.png)   

### 2.1  Quality Control using sickle  
**working directory**  
```
Genome_Assembly/
├── short_read_assembly/
       ├── 02_quality_control/
``` 
Sickle takes raw reads and outputs data with the 3’ and 5’ ends trimmed to assure that the quality of the read is high enough for assembly, it will also trim low quality reads.  

```bash
module load sickle/1.33

module load sickle/1.33

sickle pe \
	-f ../01_raw_reads/Sample_R1.fastq \
	-r ../01_raw_reads/Sample_R2.fastq \
	-t sanger \
	-o trim_Sample_R1.fastq \
	-p trim_Sample_R2.fastq \
	-s sinlges.fastq \
	-q 30 \
	-l 45 

module unload sickle/1.33
```  
The useage information on the sickle program:  
```
Usage: sickle pe [options] -f <paired-end forward fastq file> 
	-r <paired-end reverse fastq file> 
	-t <quality type> 
	-o <trimmed PE forward file> 
	-p <trimmed PE reverse file> 
	-s <trimmed singles file>    

Options:
-f, --pe-file1, Input paired-end forward fastq file
-r, --pe-file2, Input paired-end reverse fastq file
-o, --output-pe1, Output trimmed forward fastq file
-p, --output-pe2, Output trimmed reverse fastq file
-s                Singles files

Global options:
-t, --qual-type, Type of quality values
                solexa (CASAVA < 1.3)
                illumina (CASAVA 1.3 to 1.7)
                sanger (which is CASAVA >= 1.8)
-s, --output-single, Output trimmed singles fastq file
-l, --length-threshold, Threshold to keep a read based on length after trimming. Default 20
-q, --qual-threshold, Threshold for trimming based on average quality in a window. Default 20
```  
The full slurm script is called [sr_quality_control.sh](short_read_assembly/02_quality_control/sr_quality_control.sh) which can be found in the *02_quality_control/* folder.
Once you run the batch script using *sbatch* command, you will end up with the following files:   
```
02_quality_control/
├── trim_Sample_R1.fastq
├── trim_Sample_R2.fastq
└── sinlges.fastq
```  

### 2.2  Assembly  
#### 2.2a  Assembly with SOAPdenovo   
[SOAP-denovo](https://www.animalgenome.org/bioinfo/resources/manuals/SOAP.html) is a short read de novo assembler. When you do deep sequencing, and have multiple libraries, they will produce multiple sequencing files. A configuration file will let the assembler know, where to find these files. In here we will provide you with a configuration file.    

The configuration file will have global information, and then multiple library sections. For global information, right now only *max_rd_len* is included. A read which is longer than this length will be cut to this length.  

Then the library information is 


this assembler, like MaSuRCA which we will be encountering later, requires a config file to run through the data. The configuration file can be found [here]().  

```bash
#maximal read length
max_rd_len=250
[LIB]
#average insert size
avg_ins=550
#if sequence needs to be reversed
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#use only first 250 bps of each read
rd_len_cutoff=250
#in which order the reads are used while scaffolding
rank=1
# cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
# path to genes
q1=/UCHC/PublicShare/Tutorials/Assembly_Tutorial/Quality_Control/Sample_1.fastq
q2=/UCHC/PublicShare/Tutorials/Assembly_Tutorial/Quality_Control/Sample_2.fastq
q=/UCHC/PublicShare/Tutorials/Assembly_Tutorial/Quality_Control/Sample_s.fastq
```

Once the configuration file is ready, and keeping it in the working directory, you can run the SOAPdenovo using following script:  
```bash
module load SOAP-denovo/2.04

cd /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Assembly/SOAP


SOAPdenovo-63mer all -s /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Assembly/Sample.config -K 31 -R -o graph_Sample_31 1>ass31.log 2>ass31.err
SOAPdenovo-63mer all -s /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Assembly/Sample.config -K 35 -R -o graph_Sample_35 1>ass35.log 2>ass35.err
SOAPdenovo-63mer all -s /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Assembly/Sample.config -K 41 -R -o graph_Sample_41 1>ass41.log 2>ass41.err

module unload SOAP-denovo/2.04
```


