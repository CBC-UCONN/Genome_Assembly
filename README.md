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
Working directory:  
```
short_read_assembly/
├── 03_assembly/
│   ├── SOAP/
```

[SOAP-denovo](https://www.animalgenome.org/bioinfo/resources/manuals/SOAP.html) is a short read de novo assembler. When you do deep sequencing, and have multiple libraries, they will produce multiple sequencing files. A configuration file will let the assembler know, where to find these files. In here we will provide you with a configuration file.    

The configuration file will have global information, and then multiple library sections. For global information, right now only *max_rd_len* is included. A read which is longer than this length will be cut to this length.  

Followed by the global information, the library information of sequencing data should be organized under [LIB] tag.   

Each libaray section will start with [LIB] tag: following are the items which it will include.   

*   **avg_ins** : average insert size of this library or the peak value position in the insert size distribution.   
*   **reverse_seq** : This option will take value 0 or 1. It tells the assembler if the read sequences need to be complementarily reversed. Illumima GA produces two types of paired-end libraries:   
    *   forward-reverse, generated from fragmented DNA ends with typical insert size less than 500 bp;
	*   forward-forward, generated from circularizing libraries with typical insert size greater than 2 Kb;   
	The parameter "reverse_seq" should be set to indicate this:   
	0, forward-reverse; 1, forward-forward.   

*   **asm_flags=3** : This indicator decides in which part(s) the reads are used. It takes value:    
    *    1 : only contig assembly 
	*    2 : only scaffold assembly   
	*    3 : both contig and scaffold assembly   
	*    4 : only gap closure    

*    **rd_len_cutoff** : The assembler will cut the reads from the current library to this length.   
*    **rank** : it takes integer values and decides in which order the reads are used for scaffold assembly. Libraries with the same "rank" are used at the same time during scaffold assembly.    
*    **pair_num_cutoff** : This parameter is the cutoff value of pair number for a reliable connection between two contigs or pre-scaffolds.   
*    **map_len** : This takes effect in the "map" step and is the minimun alignment length between a read and a contig required for a reliable read location.   

After the above tags the reads can be added in the following fashion:  
*   It will accept two formats: FASTA or FASTQ   
	*   single end files are indicated in **f=/path-to-file/** or **q=/path-to-file/** : FASTA / FASTQ   
	*   paired reads in two FASTA sequence files are indicated by: **f1=/path-to-file/**  and **f2=/path-to-file/**   
	*   paired reads in two fastq sequences files are indicated by: **q1=/path-to-file/** and **q2=/path-to-file/**   
	*   paired reads in a single fasta sequence file is indicated by **"p="**   

following is the configuration file we are using for our run:   

```bash
#Global
#maximal read length
max_rd_len=250
#Library
[LIB]
#average insert size
avg_ins=550
#if sequence needs to be reversed in our case its; forward-reverse
reverse_seq=0
#both contig and scaffold assembly parts of the reads are used
asm_flags=3
#use only first 250 bps of each read
rd_len_cutoff=250
#in which order the reads are used while scaffolding
rank=1
# cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
#PATH to FASTQ reads 
q1=../../02_quality_control/trim_Sample_R1.fastq
q2=../../02_quality_control/trim_Sample_R2.fastq
q=../../02_quality_control/sinlges.fastq
```

The configuration file we created for this job is named as [config_file](short_read_assembly/03_assembly/config_file) and can be found in the **SOAP/** folder.   

Once the configuration file is ready, you can run the assembly job using the following command:   

```bash
module load SOAP-denovo/2.04

SOAPdenovo-63mer all \
	-s config_file \
	-K 31 \
	-p 8 \
	-R \
	-o graph_Sample_31 1>kmer31.log 2>kmer31.err 

SOAPdenovo-63mer all \
	-s config_file \
	-K 35 \
	-p 8 \
	-R \
	-o graph_Sample_35 1>kmer35.log 2>kmer35.err 

SOAPdenovo-63mer all \
	-s config_file \
	-K 41 \
	-p 8\
	-R \
	-o graph_Sample_41 1>kmer41.log 2>kmer41.err 
```

SOAPdenovo assembly options:   
```
 Usage: SOAPdenovo <command> [option]

all             do pregraph-contig-map-scaff in turn

  -s <string>    configFile: the config file
  
  -K <int>       kmer(min 13, max 63): kmer size
  -p <int>       n_cpu: number of cpu for use
  -R (optional)  resolve repeats by reads
  -o <string>    outputGraph: prefix of output graph file name     
 ```

> NOTE
> A k-mer is a set of nucleotides, k is the number of nucleotides in that set. It is a crucial parameter in most de Brujin Graph assemblers and assemblers work with the highest accuracy if the k-mer size estimation is accurate.

The above script is called [SOAPdenovo.sh](/short_read_assembly/03_assembly/SOAPdenovo.sh) and can be found in the SOAP/ directory. The script can be run using the `sbatch` command.  

It will produce bunch of files, and we are interested in the each k-mer scafold sequences produced, at the end of each run. These are the files which we will be used to asses the quality of our assembly.  


```
SOAP
├── graph_Sample_31.scafSeq
├── graph_Sample_35.scafSeq
└── graph_Sample_41.scafSeq
```   


#### 2.2b  Assembly with SPAdes   
Workding directory:
```
short_read_assembly/
└── 03_assembly/
    └── SPAdes/  
```

[SPAdes](http://cab.spbu.ru/software/spades/) - St. Petersburg genome assembler is a toolkit containing assembly pipelines. When SPAdes was initally designed it was for for small genonmes, like bacterial, fungal and other small genomes. SPAdes is not intened for larger genomes.SPAdes has different pipe-lines present and if you want to check them out please visit the [SPAdes web site](http://cab.spbu.ru/software/spades/) or its [git-hub site](https://github.com/ablab/spades/blob/spades_3.14.0/README.md).  

Instead of manually selecting k-mers, SPAdes automatically selects k-mers based off the maximum read length data of your input. This is a called a de Bruijn graph based assembler, meaning that it assigns (k-1)-mers to nodes and every possible matching prefix and suffix of these nodes are connected with a line.  


SPAdes takes as input paired-end reads, mate-pairs and single (unpaired) reads in FASTA and FASTQ. In here we are using paired-end reads and the left and right reads are held in two files and they will be taken into the assembler in the same order as in the files.  

Command line options we used:   
```bash
module load SPAdes/3.13.0

spades.py \
	-1 ../../02_quality_control/trim_Sample_R1.fastq \
	-2 ../../02_quality_control/trim_Sample_R2.fastq \
	-s ../../02_quality_control/sinlges.fastq \
	--careful \
	--threads 8 \
	--memory 30 \
	-o .
```

Basic SPAdes command line would look like:  
`spades.py [options] -o <output_dir>`   

where the options we used:
```
Input
-1	file with forward paired-end reads
-2	file with reverse paired-end reads
-s	file with unpaired reads

--careful    tries to reduce number of mismatches and short indels
--threads	 number of threads
--memory	 RAM limit for SPAdes in Gb (terminates if exceeded) defaul is 250
```  
> **NOTE**  
> its very important to make sure you match the number of theads and the memory asked in the options section is matched with the SLURM header part of your script.  

The full script for our SPAdes run is called [SPAdes.sh](/short_read_assembly/03_assembly/SPAdes/SPAdes.sh) and it can be found in the **03_assembly/SPAdes/** directory.  

Once the assembly script is ran, it will produce bunch of files together with the final scafold file which is called, scaffolds.fasta. We will be using this final scaffold file to analyze the SPAdes assembly run. 
```
SPAdes/
├── scaffolds.fasta
 
```


#### 2.2c Assembly with MaSuRCA   

Working directory will be:  
```
03_assembly/
├── MaSuRCA/
```
[MaSuRCA](https://github.com/alekseyzimin/masurca) (**Ma**ryland **Su**per **R**ead **C**abog **A**ssembler) is a combination of a De Bruijn graph and an Overlap-Layout-Consensus model. The Overlap-Layout-Consensus model consists of three steps, Overlap, which is the process of overlapping matching sequences in the data, this forms a long branched line. Layout, which is the process of picking the least branched line in from the overlap sequence created earlier, the final product here is called a contig. Consensus is the process of lining up all the contigs and picking out the most similar nucleotide line up in this set of sequences (OIRC).   

When running [MaSuRCA](http://www.genome.umd.edu/masurca.html), there are few things you should keep in mind. This assembler, **DOES NOT** require a preprocessing step, such as trimming, cleaning or error correction step; you will directly feed the raw reads.   

The first step is to create a configuration file. A sample file can be copied from the MaSuRCA instalation directory. The following command will copy it to your working folder.  
```bash
cp $MASURCA/sr_config_example.txt config_file
```  
The configuration file, will contain the location of the compiled assembler, the location of the data and some parameters. In most cases you only need to change the path to your read files.   

Second step is to run the masurca script which will create a shell script called `assembly.sh` using the configuration file.  

Then the final step is to run this `assembly.sh` script, which will create the scaffolds.   

Lets look at the configuration file, which contain two sections, DATA and PARAMERTERS, and each section concludes with END section.  

In the DATA section:
```bash
DATA
#Illumina paired end reads supplied as <two-character prefix> <fragment mean> <fragment stdev> <forward_reads> <reverse_reads> 
PE= pe 180 20  ../../01_raw_reads/Sample_R1.fastq ../../01_raw_reads/Sample_R2.fastq 
END
```   
DATA section is where you should specify the input data for the assembler. Each library line should with the appropiate read type, `eg: PE, JUMP, OTHER`.  In the above DATA section we have specified Illumina paired end reads.  
`PE = two_letter_prefix mean stdev /path-to-forward-read /path-to-reverse-read`   

 `mean` = is the library insert average length  
 `stdev` = is the stanard deviation. It this is not known set it as 15% of the `mean`. If the reverse read is not avaliable do not specify this.      
If you are interested in other types of reads and how to include them in the DATA section, more information can be found in the [MaSuRCA git page](https://github.com/alekseyzimin/masurca).  


In the PARAMETERS section:  
```bash
PARAMETERS
#PLEASE READ all comments to essential parameters below, and set the parameters according to your project
#set this to 1 if your Illumina jumping library reads are shorter than 100bp
EXTEND_JUMP_READS=0
#this is k-mer size for deBruijn graph values between 25 and 127 are supported, auto will compute the optimal size based on the read data and GC content
GRAPH_KMER_SIZE = auto
#set this to 1 for all Illumina-only assemblies
#set this to 0 if you have more than 15x coverage by long reads (Pacbio or Nanopore) or any other long reads/mate pairs (Illumina MP, Sanger, 454, etc)
USE_LINKING_MATES = 0
#specifies whether to run the assembly on the grid
USE_GRID=0
#specifies grid engine to use SGE or SLURM
GRID_ENGINE=SGE
#specifies queue (for SGE) or partition (for SLURM) to use when running on the grid MANDATORY
GRID_QUEUE=all.q
#batch size in the amount of long read sequence for each batch on the grid
GRID_BATCH_SIZE=500000000
#use at most this much coverage by the longest Pacbio or Nanopore reads, discard the rest of the reads
#can increase this to 30 or 35 if your reads are short (N50<7000bp)
LHE_COVERAGE=25
#set to 0 (default) to do two passes of mega-reads for slower, but higher quality assembly, otherwise set to 1
MEGA_READS_ONE_PASS=0
#this parameter is useful if you have too many Illumina jumping library mates. Typically set it to 60 for bacteria and 300 for the other organisms 
LIMIT_JUMP_COVERAGE = 300
#these are the additional parameters to Celera Assembler.  do not worry about performance, number or processors or batch sizes -- these are computed automatically. 
#CABOG ASSEMBLY ONLY: set cgwErrorRate=0.25 for bacteria and 0.1<=cgwErrorRate<=0.15 for other organisms.
CA_PARAMETERS =  cgwErrorRate=0.15
#CABOG ASSEMBLY ONLY: whether to attempt to close gaps in scaffolds with Illumina  or long read data
CLOSE_GAPS=1
#auto-detected number of cpus to use, set this to the number of CPUs/threads per node you will be using
NUM_THREADS = 16
#this is mandatory jellyfish hash size -- a safe value is estimated_genome_size*20
JF_SIZE = 200000000
#ILLUMINA ONLY. Set this to 1 to use SOAPdenovo contigging/scaffolding module.  Assembly will be worse but will run faster. Useful for very large (>=8Gbp) genomes from Illumina-only data
SOAP_ASSEMBLY=0
#Hybrid Illumina paired end + Nanopore/PacBio assembly ONLY.  Set this to 1 to use Flye assembler for final assembly of corrected mega-reads.  A lot faster than CABOG, at the expense of some contiguity. Works well even when MEGA_READS_ONE_PASS is set to 1.  DO NOT use if you have less than 15x coverage by long reads.
FLYE_ASSEMBLY=0
END
```

Once the configuration file is set up you can run MaSuRCA using:
```bash
module load MaSuRCA/3.3.4 

masurca config_file

bash assemble.sh

module unload MaSuRCA/3.3.4
```

The full script for running MaSuRCA is called [MASuRCA.sh](short_read_assembly/03_assembly/MaSuRCA/MASuRCA.sh) and can be found in the 03_assembly/MaSuRCA/ folder.   

Final assembly scaffolds can be found under the **CA/** folder:  
```
MaSuRCA/
├── CA/
│   ├── final.genome.scf.fasta
```

### 2.3  Quality Assesment  

Overview: Once you have the contigs/scaffolds assembled, next step is to check the quality of these. In here we will use 4 methods to check, which will be depicted in the following:  

![](/images/quality_assesment_overview-1.png)  

