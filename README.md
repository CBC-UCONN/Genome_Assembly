# Genome Assembly   

This repository is a usable, publicly available tutorial. All steps have been provided for the UConn CBC Xanadu cluster here with appropriate headers for the Slurm scheduler that can be modified simply to run.  Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should use sbatch scriptname after modifying the script for each stage.  Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs.  If you are new to Linux, please use [this](https://bioinformatics.uconn.edu/unix-basics) handy guide for the operating system commands.  In this guide, you will be working with common bioinformatic file formats, such as [FASTA](https://en.wikipedia.org/wiki/FASTA_format), [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), [SAM/BAM](https://en.wikipedia.org/wiki/SAM_(file_format)), and [GFF3/GTF](https://en.wikipedia.org/wiki/General_feature_format). You can learn even more about each file format [here](https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/). If you do not have a Xanadu account and are an affiliate of UConn/UCHC, please apply for one **[here](https://bioinformatics.uconn.edu/contact-us/)**.   

**Contents**    

1.   [Overview](#1-overview)  
2.   [Short Read Genome Assembly](#2-short-read-genome-assembly)
      1.  [Quality Control](#21--quality-control-using-sickle)  
      2.  [Assembly](#22--assembly)   
      3.  [Quality Assesment](#23--quality-assesment)  
3.   [Long Read Genome Assembly](#3-long-read-genome-assembly)  
4.   [Hybrid Genome Assembly](#4-hybrid-genome-assembly)  


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

 In this study we are using short reads produced by Illumina MiSeq instrument. This data refers to a paired end reads with a read length of 250bp. The Nextera preperation method is sensitive and can produce lots of small fragments, shorter than 500bp. This resutls in R1 and R2 actually overlaping each other.    

 ![](images/Nextera_library.jpg)  

This has a estimated genome size is 3Mb and a inital coverage of 84**X**. As a problem you can try to calculate the coverage of this data set. Check out on how to calculate the coverage [here](coverage.md).  


### 2.1  Quality Control  


#### Quality Check of Reads using FASTQC  
In here we will use the FASTQC package to check the quality of the reads before and after trimming.   

Quality check of raw reads:  
```
mkdir -p RAWfastqc_OUT
fastqc -o ./RAWfastqc_OUT ../01_raw_reads/Sample_R1.fastq ../01_raw_reads/Sample_R2.fastq
```  

Quality check of the reads show a adapter contamination of the reads from Nextera sequences, which need to be removed.  

#### Trimming of reads using Trimmomatic  

**working directory**  
```
Genome_Assembly/
├── short_read_assembly/
       ├── 02_quality_control/
```  

Trimmomatic is used to remove low quality and adapter sequences. The program can be used on single end as well as paired end sequences. Since we have identified the raw reads have Nextera sequences present we have the option to give specific sequences for the program to trim using the `ILLUMINACLIP` command. More information on the specific options can be found in [Trimmomatic website](http://www.usadellab.org/cms/?page=trimmomatic). 

```
java -jar $Trimmomatic PE -threads 4 \
        ../01_raw_reads/Sample_R1.fastq \
        ../01_raw_reads/Sample_R2.fastq \
        trim_Sample_R1.fastq trim_Sample_R1_singles.fastq \
        trim_Sample_R2.fastq trim_Sample_R2_singles.fastq \
        ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 \
        SLIDINGWINDOW:4:25 MINLEN:45
```  

The complete slurm scrip is called [sr_quality_control.sh](short_read_assembly/02_quality_control/sr_quality_control.sh). It will produce four fastq files, where two will contain trimmed reads of forward and reverse sequences and two files of single fastq sequences where only one sequence (forward or reverse) pass the filtering criteria.   

```
02_quality_control/
├── trim_Sample_R1.fastq
├── trim_Sample_R1_singles.fastq
├── trim_Sample_R2.fastq
└── trim_Sample_R2_singles.fastq
```  

fter trimming the coverage will be ~80**X** 


#### Quality check of trimmed reads:  
```
mkdir -p TRIMfastqc_OUT
fastqc -o ./TRIMfastqc_OUT ./trim_Sample_R1.fastq ./trim_Sample_R2.fastq
```  

FASTQC produces a HTML file with stats about your reads. You can download these HTML files to your local computer to view them using the `transfer.cam.uchc.edu` submit node, which facilitate file transfer.   




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

The configuration file we created for this job is named as [config_file](short_read_assembly/03_assembly/SOAP/config_file) and can be found in the **SOAP/** folder.   

Once the configuration file is ready, you can run the assembly job using the following command:   


```
SOAPdenovo-127mer all \
        -s config_file \
        -K 31 \
        -p 8 \
        -R \
        -o graph_Sample_31 1>ass31.log 2>ass31.err 

SOAPdenovo-127mer all \
        -s config_file \
        -K 71 \
        -p 8 \
        -R \
        -o graph_Sample_71 1>ass71.log 2>ass71.err 

SOAPdenovo-127mer all \
        -s config_file \
        -K 101 \
        -p 8\
        -R \
        -o graph_Sample_101 1>ass101.log 2>ass101.err 
```

SOAPdenovo assembly options:   
```
 Usage: SOAPdenovo-127mer <command> [option]

all             do pregraph-contig-map-scaff in turn

  -s <string>    configFile: the config file
  
  -K <int>       kmer(min 13, max 127): kmer size
  -p <int>       n_cpu: number of cpu for use
  -R (optional)  resolve repeats by reads
  -o <string>    outputGraph: prefix of output graph file name     
 ```

> NOTE
> A k-mer is a set of nucleotides, k is the number of nucleotides in that set. It is a crucial parameter in most de Brujin Graph assemblers and assemblers work with the highest accuracy if the k-mer size estimation is accurate.

The above script is called [SOAPdenovo.sh](short_read_assembly/03_assembly/SOAP/SOAPdenovo.sh) and can be found in the SOAP/ directory. The script can be run using the `sbatch` command.  

It will produce bunch of files, and we are interested in the each k-mer scafold sequences produced, at the end of each run. These are the files which we will be used to asses the quality of our assembly.  


```
SOAP
├── graph_Sample_31.scafSeq
├── graph_Sample_71.scafSeq
└── graph_Sample_101.scafSeq
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
PE= pe 480 20  ../../01_raw_reads/Sample_R1.fastq ../../01_raw_reads/Sample_R2.fastq 
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
CA_PARAMETERS =  cgwErrorRate=0.25
#CABOG ASSEMBLY ONLY: whether to attempt to close gaps in scaffolds with Illumina  or long read data
CLOSE_GAPS=1
#auto-detected number of cpus to use, set this to the number of CPUs/threads per node you will be using
NUM_THREADS = 8
#this is mandatory jellyfish hash size -- a safe value is estimated_genome_size*20
JF_SIZE = 300000000
#ILLUMINA ONLY. Set this to 1 to use SOAPdenovo contigging/scaffolding module.  Assembly will be worse but will run faster. Useful for very large (>=8Gbp) genomes from Illumina-only data
SOAP_ASSEMBLY=0
#Hybrid Illumina paired end + Nanopore/PacBio assembly ONLY.  Set this to 1 to use Flye assembler for final assembly of corrected mega-reads.  A lot faster than CABOG, at the expense of some contiguity. Works well even when MEGA_READS_ONE_PASS is set to 1.  DO NOT use if you have less than 15x coverage by long reads.
FLYE_ASSEMBLY=0
END
```

The full configuration file in our run is called [config_file](short_read_assembly/03_assembly/MaSuRCA/config_file) and it can be found in the MaSuRCA/ directory. Once the configuration file is set up you can run MaSuRCA using:
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

## 2.3  Quality Assesment  

Once we have the contigs/scaffolds assembled, next step is to check the quality of these. In here we will use 4 methods to check, which will be depicted in the following:  

![](/images/short_read_qc.jpg)  


### 2.3.a Assembly Statistics with QUAST   

[QUAST](http://quast.sourceforge.net/quast) will be used to evaluate genome assemblies. We will be using the program QUAST which will give us the number of contigs, total length and N50 value; the data we are most interested in. A good assembly would have small number of contigs, a total length that makes sense for the specific species, and a large N50 value. N50 is a measure to describe the quality of assembled genomes fragmented in contigs of different length. The N50 is the minimum contig length needed to cover 50% of the genome.  

Working directory:  
```
short_read_assembly/
└── 04_quast
```



*   For SOAPdenovo scaffolds we will be using:  
```
module load quast/5.0.2

quast.py ../03_assembly/SOAP/graph_Sample_31.scafSeq \
        --threads 8 \
        -o SOAP_31

quast.py ../03_assembly/SOAP/graph_Sample_71.scafSeq \
        --threads 8 \
        -o SOAP_71

quast.py ../03_assembly/SOAP/graph_Sample_101.scafSeq \
        --threads 8 \
        -o SOAP_101
```
The full script is called [quast_SOAP.sh](short_read_assembly/04_quast/quast_SOAP.sh) and can be found in the 04_quast directory.  

*  For SPAdes; the command we use would be like:  
```
module load quast/5.0.2

quast.py ../03_assembly/SPAdes/scaffolds.fasta \
	--threads 8 \
	-o SPAdes
```  
The full script is called [quast_SPAdes.sh](short_read_assembly/04_quast/quast_SPAdes.sh) and can be found in the 04_quast directory.  

*   For MaSuRCA; the command we use would be like:  
```
module load quast/5.0.2

quast.py ../03_assembly/MaSuRCA/CA/final.genome.scf.fasta \
	--threads 8 \
	-o MaSuRCA
```   

The full script is called [quast_MaSuRCA.sh](short_read_assembly/04_quast/quast_MaSuRCA.sh), and can be found in the 04_quast directory.   

General command for executing quast would be like:  
`quast.py [options] <files_with_contigs>`    

Options:  
```
-o  --output-dir  Directory to store all result files
-t  --threads     Maximum number of threads [default: 25% of CPUs]

```  

These are the minimum command that we used, there are many other options that you can use with the software and if you are interested please have a look in to the [QUAST manual](http://quast.sourceforge.net/docs/manual.html).

Once executed these scripts using the `sbatch` command, you will end of with basic evaluation of the assemblies. These statics can be found in each new folder you created:
```
04_quast/
├── MaSuRCA
│   └── report.txt
├── SOAP_31
│   └── report.txt
├── SOAP_71
│   └── report.txt
├── SOAP_101
│   └── report.txt
└── SPAdes
    └── report.txt
```

 
|             |  SOAP-31     |  SOAP-71     |  SOAP-101     |  SPAdes   |   MaSuRCA      |    
 ------------ |:---------: | :---------: | ---------: | --------- | ---------- |  
 contigs (>= 0 bp)    | 1507  | 4252  |893   |101   | 86 | 
 contigs (>= 1000 bp) | 249   |  159  | 154   | 52   | 77 | 
 contigs (>= 5000 bp) | 158   |  113  | 106   | 41   | 69 | 
 contigs (>= 10000 bp) | 115  |  90  | 85   | 35   | 62 | 
 contigs (>= 25000 bp) | 47   |  41   | 39    | 27   | 40 | 
 contigs (>= 50000 bp) | 5    |  14   | 11    | 19   | 18 | 
 Total length (>= 0 bp)| 3743924 |  3656745 | 3008494 |  2885853 | 2825384 | 
 Total length (>= 1000 bp)| 3554783 | 3039213 | 2861233 | 2870111 | 2817661 | 
 Total length (>= 5000 bp)| 3302110 | 2912904 | 2727264 | 2843776 | 2796700 | 
 Total length (>= 10000 bp)| 2985682 | 2732622 | 2578612 | 2799884 | 2745806 | 
 Total length (>= 25000 bp)| 1871285 |  1947727 | 1828836 | 2664507 | 2388905 | 
 Total length (>= 50000 bp)| 371124 |  982749 | 814206 | 2373293 | 1610768 | 
 **no. of contigs**     |  276  | 173 | 169 | 60  | 86 | 
 Largest contig     | 103125 |  98502 | 173926 | 255651 | 220219 |
 Total length       | 3574101 |  3049168 | 2871343 | 2875218 | 2825384 | 
 GC (%)      | 32.44  |  32.55 | 32.62 | 32.65 | 32.75 | 
 **N50**         | 26176  |  37961 | 36337 | 149694 | 54020 |
 N75         | 14642  |  17093 | 18263 | 61620 | 32399 | 
 L50         |  44  |  27 | 26 | 8 | 15 | 
 L75         |  91  |  57 | 54 | 15 | 31 |  



According to our requirements regarding n50 and contigs it would appear that the best assembly perfromed was via SPAdes. (N50 value indicates that half the genome is assembled on contigs/scaffolds of length N50 or longer)


### 2.3.b Read Alignment with Bowtie2   

[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-aligner) is a tool you would use for comparitive genomics via alignment. Bowtie2 takes a Bowtie 2 index and a set of sequencing read files and outputs a set of alignments in a SAM file. Alignment is the process where we discover how and where the read sequences are similar to a reference sequence. An alignment is a way of lining up the characters in the read with some characters from the reference to reveal how they are similar.  

Alignemnt is a method of doing a educated guess, as where the read originated with respect to the reference genome. It is not always possibe to determine this with clarity.   

Bowtie2 in our case takes read sequences and aligns them with long reference sequences. Since this is de novo assembly you will take the data from the assemblies you have and align them back to the raw read data. You want to use unpaired data.  

Our working directory will be:  
```
short_read_assembly/
├── 05_bowtie2
```  

First step in aligning reads with bowtie2 is to make a index of the genome we assembled. This can be done with the `bowtie2-build` command.  

```bash
bowtie2-build [options] <reference-index> <index-base-name>  
```  

Once you build the index, next step would be to align the reads to the genome using the index you build.   

```bash
bowtie2 [options] -x <bt-index> -1 <paired-reads> -2 <paired-reads> -S <SAM-output> \
	--threads 8 2>output.err
```   

In here we will direct the error file output which will contain the alignment statistics.   

Since we are using three de novo methods to construct genomes, we will now try to see how each of them are aligning its reads to the genome constructed.   

*   SOAP   
```bash
mkdir -p SOAP_31_index SOAP_71_index SOAP_101_index

module load bowtie2/2.3.5.1

## SOAP_31
bowtie2-build \
        --threads 8 \
        ../03_assembly/SOAP/graph_Sample_31.scafSeq SOAP_31_index/SOAP_31_index

bowtie2 -x SOAP_31_index/SOAP_31_index \
        -1 ../01_raw_reads/Sample_R1.fastq -2 ../01_raw_reads/Sample_R2.fastq \
        -S SOAP_31.bowtie2.sam \
        --threads 8 2>SOAP_31.err

## SOAP_71
bowtie2-build \
        --threads 8 \
        ../03_assembly/SOAP/graph_Sample_71.scafSeq SOAP_71_index/SOAP_71_index

bowtie2 -x SOAP_71_index/SOAP_71_index \
        -1 ../01_raw_reads/Sample_R1.fastq -2 ../01_raw_reads/Sample_R2.fastq \
        -S SOAP_71.bowtie2.sam \
        --threads 8 2>SOAP_71.err

## SOAP_101 
bowtie2-build \
        --threads 8 \
        ../03_assembly/SOAP/graph_Sample_101.scafSeq SOAP_101_index/SOAP_101_index

bowtie2 -x SOAP_101_index/SOAP_101_index \
        -1 ../01_raw_reads/Sample_R1.fastq -2 ../01_raw_reads/Sample_R2.fastq \
        -S SOAP_101.bowtie2.sam \
        --threads 8 2>SOAP_101.err

```  

The full slurm script is called [bowtie2_SOAP.sh](/short_read_assembly/05_bowtie2/bowtie2_SOAP.sh), and can be found in the 05_bowtie2 directory.  


*   SPAdes  
```bash
# SPAdes
mkdir -p SPAdes_index

module load bowtie2/2.3.5.1

bowtie2-build \
        --threads 8 \
        ../03_assembly/SPAdes/scaffolds.fasta SPAdes_index/SPAdes_index

bowtie2 -x SPAdes_index/SPAdes_index \
        -1 ../01_raw_reads/Sample_R1.fastq -2 ../01_raw_reads/Sample_R2.fastq \
        -S SPAdes.bowtie2.sam \
        --threads 8 2>SPAdes.err
```

The full slurm script for running bowtie2 for genome created with  is called [bowtie2_SPAdes.sh](/short_read_assembly/05_bowtie2/bowtie2_SPAdes.sh), which can be found in **05_bowtie2/** directory.  

*   MaSuRCA   
```bash
mkdir -p MaSuRCA_index

module load bowtie2/2.3.5.1

bowtie2-build \
        --threads 8 \
        ../03_assembly/MaSuRCA/CA/final.genome.scf.fasta MaSuRCA_index/MaSuRCA_index

bowtie2 -x MaSuRCA_index/MaSuRCA_index \
        -1 ../01_raw_reads/Sample_R1.fastq -2 ../01_raw_reads/Sample_R2.fastq \
        -S MaSuRCA.bowtie2.sam \
        --threads 8 2>MaSuRCA.err
```  

The full slurm script for running bowtie2 for genome created with MaSuRCA is called [bowtie2_MaSuRCA.sh](/short_read_assembly/05_bowtie2/bowtie2_MaSuRCA.sh), which can be found in 05_bowtie2/ directory.   

As shown for the SOAP_31: it will create the following files and folders for each above (SOAP_35, SOAP_45, SPAdes, MaSuRCA) runs:  
```
05_bowtie2/
├── SOAP_31_index/
│   ├── SOAP_31_index.1.bt2
│   ├── SOAP_31_index.2.bt2
│   ├── SOAP_31_index.3.bt2
│   ├── SOAP_31_index.4.bt2
│   ├── SOAP_31_index.rev.1.bt2
│   └── SOAP_31_index.rev.2.bt2
├── SOAP_31.bowtie2.sam
└── SOAP_31.err
``` 
We have only shown the output directory and files for SOAP_131 case, as an example.   

The alingment data will be in `*.err` files associated with each run.  
```
05_bowtie2/
├── SOAP_31.err
├── SOAP_71.err
├── SOAP_101.err
├── SPAdes.err
└── MaSuRCA.err
```    

Following table summerizes the alignment results with bowtie2.  


|             |  SOAP-31     |  SOAP-71     |  SOAP-101     |  SPAdes   |   MaSuRCA      |    
 ------------ |:---------: | :---------: | ---------: | --------- | ---------- |   
reads   |  504951  |   504951   |  504951  |  504951  |  504951  |  
paired   |  504951   |   504951   |  504951  |  504951  |  504951  |  
aligned concordantly 0 times  |  76.96%   |   59.41%  |  45.56% |  42.00%  |  43.80%  |  
aligned concordantly exactly 1 time   |  23.04%   |  40.51%  |   54.09%  |  57.32%  |  53.84%  |   
aligned concordantly >1 times   |   0.00%   |   0.09%  | 0.35% |  0.68%  |  2.36%  |  
overall alignment rate   |  44.66%   |   66.63%  |  82.63% |  86.77%  |   83.89%  |  


Between de-novo assemblies it shows that, assembly done with SPAdes have a good overall alignment rate, and higher number of reads would match exactly 1 time to the reference genome.    


### 2.3.c BUSCO evaluation: Assessing Genome Assembly and Annotation Completeness    

Here, we describe the use of the [BUSCO](https://busco.ezlab.org/) tool suite to assess the completeness of genomes, gene sets, and transcriptomes, using their gene content as a complementary method to common technical metrics.       

To access the quality and completeness of a assembly, different matrices can be used. Such as alignment ratios and contig/scaffold length distributions as N50 will reflect the assembly completeness. The above technics ignore the biological aspects regarding the gene content. This would be an important aspect where you have transcriptomic data. You can test comprehensiveness of the gene set by aligning the transcripts to the assembly to assess how much is aligned. However, aligning the spliced transcripts to genomic regions can be problematic and it depends on the tools and parameters used for mapping. This leads us to look for another alternative, and [BUSCO](https://busco.ezlab.org/): **B**enchmarking **U**niversal **S**ingle-**C**opy **O**rthologs provides a way of assessing genome assemblies.   

BUSCO data set is made up of protein coding genes from the OrthoDB orthologus groups, that contain single copy genes in that group. The consensus sequences profiles are built from multiple alignment across the species group. As more species sequenced the BUSCO data set will be updated with more lineages. This program assists with checking assemblies, gene sets, annotations, and transcriptomes to see if they appear complete, using their gene content as a complementary method to common technical metrics. It does this by taking an orthologous gene set of your species of interest and comparing it back to the genome of interest, taking into consideration possibly evolutionary changes.    

We will now try to evaluate our three assembled genomes using BUSCO. For this we will be using the BUSCO database which is downloaded in Xanadu cluster (`/isg/shared/databases/BUSCO/`). When you are using the BUSCO database make sure you are using the latest database by checking the BUSCO database page. In this tutorial we are using the current database version which is ‘**odb10**’, which is compatible with BUSCO software version 4.   

Our working directory will be:    
```
short_read_assembly/
└── 06_busco
```  

Following will be the commands which will be used for evaluating SOAP, SPAdes, MaSuRCA assemblies using BUSCO.   

*   SOAP:  

```bash  
busco -i ../03_assembly/SOAP/graph_Sample_31.scafSeq \
        -o SOAP_31 -l /isg/shared/databases/BUSCO/odb10/bacteria_odb10 -m genome


busco -i ../03_assembly/SOAP/graph_Sample_71.scafSeq \
        -o SOAP_71 -l /isg/shared/databases/BUSCO/odb10/bacteria_odb10 -m genome


busco -i ../03_assembly/SOAP/graph_Sample_101.scafSeq \
        -o SOAP_101 -l /isg/shared/databases/BUSCO/odb10/bacteria_odb10 -m genome
```   

The full script for evaluating the three SOAP assemblies is called, [busco_SOAP.sh](short_read_assembly/06_busco/busco_SOAP.sh) which can be found in `short_read_assembly/06_busco` folder.   


*   SPAdes:  

```bash
busco -i ../03_assembly/SPAdes/scaffolds.fasta \
	-o SPAdes \
	-l /isg/shared/databases/BUSCO/odb10/bacteria_odb10 \
	-m genome
```   

The full script for evaluating the SPAdes assembly is called, [busco_SPAdes.sh](short_read_assembly/06_busco/busco_SPAdes.sh) which can be found in `short_read_assembly/06_busco` folder.  

*   MaSuRCA:  
```bash
busco -i ../03_assembly/MaSuRCA/CA/final.genome.scf.fasta \
	-o MaSuRCA \
	-l /isg/shared/databases/BUSCO/odb10/bacteria_odb10 \
	-m genome
```   

The full script for evaluating the MaSuRCA assembly is called, [busco_MaSuRCA.sh](short_read_assembly/06_busco/busco_MaSuRCA.sh) which can be found in `short_read_assembly/06_busco` folder.   

So the general command for evaluating the assembly will be:  
```
usage: busco -i [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] -m [MODE] [OTHER OPTIONS]  
```  

The command options which we will be using:  
```
-i FASTA FILE   Input sequence file in FASTA format
-l LINEAGE      Specify the name of the BUSCO lineage
-o OUTPUT       Output folders and files will be labelled with this name
-m MODE         BUSCO analysis mode
					- geno or genome, for genome assemblies (DNA)
					- tran or transcriptome, for transcriptome assemblies (DNA)
					- prot or proteins, for annotated gene sets (protein)
```  

In our run we will be evaluating a genome assembly, which is nucleotide sequence in contigs or scaffolds. In the genome run we will be using the options specified above. In this example we are only selecting to scan over the bacterial database in BUSCO. Make sure to use the database which suites your assembly, to check the other databases in the BUSCO database please check the BUSCO website or to check the current databases simply use the following command after loading the module: `busco --list-datasets`. In Xanadu cluster, we have downloaded the BUSCO databases so you do not have to download by your self, to check the path associated to calling the lineages please check the [Xanadu database page](https://bioinformatics.uconn.edu/databases/).    

Once you have executed the above commands the BUSCO output will contain a printed score to the standard output file and also a comprehensive output to the ‘short_summary.specific.bacteria_odb10.[OUTPUT_folder_name].txt’ file inside the output destination specified in your command. In our case:   

```
06_busco/
├── SOAP_31/
│   └── short_summary.specific.bacteria_odb10.SOAP_31.txt
├── SOAP_71/
│   └── short_summary.specific.bacteria_odb10.SOAP_71.txt
├── SOAP_101/
│   └── short_summary.specific.bacteria_odb10.SOAP_101.txt
├── SPAdes/
│   └── short_summary.specific.bacteria_odb10.SPAdes.txt
└── MaSuRCA/
    └── short_summary.specific.bacteria_odb10.MaSuRCA.txt
```   

The the above output will contain:  
```
SOAP_31
C:45.2%[S:45.2%,D:0.0%],F:39.5%,M:15.3%,n:124

SOAP_71
C:68.5%[S:68.5%,D:0.0%],F:26.6%,M:4.9%,n:124

SOAP_101
C:96.8%[S:96.8%,D:0.0%],F:3.2%,M:0.0%,n:124

SPAdes
C:100.0%[S:100.0%,D:0.0%],F:0.0%,M:0.0%,n:124

MaSuRCA
C:100.0%[S:100.0%,D:0.0%],F:0.0%,M:0.0%,n:124 
```
These BUSCO output will produce its output using a scoring scheme: 
**C**:complete [**S**:single-copy, **D**:duplicated], **F**:fragmented, and **M**:missing and the total BUSCO genes are indicated in **n:**.   

![](images/busco_assesment_results2.png)  

To judge the score, you need to consider the type of sequence first. A model organism with a reference genome often will reach a score of 95% or above  as a complete score and a non-model organisms can reach a score from 50% to 95% complete. This alone will not give an idea on how good the assembly is, as you need to look at the assembly and the annotation results together to make a judgement.   

In *full_table.tsv* file it will contain the detailed list of BUSCO genes and their predicted status in the genome. These files can be found in:   
```
06_busco/
├── SOAP_31
│   └── run_bacteria_odb10
│       └── full_table.tsv
├── SOAP_71
│   └── run_bacteria_odb10
│       └── full_table.tsv
├── SOAP_101
│   └── run_bacteria_odb10
│       └── full_table.tsv
├── SPAdes
│   └── run_bacteria_odb10
│       └── full_table.tsv
└── MaSuRCA
    └── run_bacteria_odb10
        └── full_table.tsv
```  

If you look into these files you will find the above information. As an example we are showing few lines of the output of SPAdes *full_table.tsv* file.

```
# Busco id      Status         Sequence                          Gene Start      Gene End        Score   Length
4421at2         Complete        NODE_19_length_52925_cov_29.157601_19   21444   25067   1718.7  1051    https://www.orthodb.org/v10?query=4421at2       DNA-directed RNA polyme
9601at2         Complete        NODE_19_length_52925_cov_29.157601_20   25204   28755   1299.7  812     https://www.orthodb.org/v10?query=9601at2       DNA-directed RNA polyme
26038at2        Complete        NODE_5_length_175403_cov_22.338455_89   93427   95616   332.5   404     https://www.orthodb.org/v10?query=26038at2      phosphoribosylf
 
```  

Other than the options we used there are many other options you can use depending on your data. More detailed view of the options are discribed in the BUSCO manual as well as in the paper: [PMID: 31020564](https://www.ncbi.nlm.nih.gov/pubmed/31020564).  



## 3. Long Read Genome Assembly   


![](images/longread_pipeline.jpg)


## 3.1 Introduction 
Long read sequencing in general has the potential to overcome the limitations that short read sequencing methods process. Long read sequencing has few main advantages over the short read technologies today.  Notably, these long sequences are produced from single DNA molecules, and the sequencing happens in real time and both sequencing and library preparation does not need the PCR amplification, which reduces the PCR based errors. This reduce the downstream errors occurring from copy errors, sequences dependent biases. Since this keeps the DNA in its native state, this gives the opportunity to detect DNA modifications such as methylation modifications. Two of the most prominent long read sequencing methods are developed by [Pacific Biosciences (PacBio)](https://www.pacb.com/) and [Oxford Nanopore Technologies (ONT)](https://nanoporetech.com/). Both single  molecule real-time (SMRT) sequencing by PacBio and nanopore sequencing by ONT, use single molecule sequencing techniques.   

In this tutorial we will focus on long reads produced by the **nanopore sequencer** (**ONT**). The basic concept of the method was developed, on the basis of understanding the change of potential patterns, when ions passes though a channel with a potential gradient. The method was developed, to identify the signal when a single DNA molecule is passed through a very narrow channel. In a sequencing run, when a DNA molecule is passed through the channel, the signal generated reflect the modulation of the ionic current in the pore. The signals produced, are stored in ‘FAST5’ format, a specialization of the HDF5 format.   

During our workflow, we will take basecalled reads and assemble the long reads using three assemblers (Flye, Shasta). Then purge Haplotags will be used to assure that the contigs that are assembled are not being combined with the Haplotig of that sequence. After this, the assembly will be polished via Nanopolish. The quality of the genome assembled, will be assessed using QUAST.  





### 3.2 Base Calling (Guppy)    

Raw reads produced by the Oxford nanopore sequencer are stored in fast5 format. Sequencing information is stored in these files as electric signal, needs to be converted into nucleic acid (A, T, G and C) format, which we can understand. This process is done by a basecalling software, called Guppy and this process is done my the sequencing center once the sequencing is done. Once the base calling is done the reads will be filtered into two groups as passed and failed. Then the passed reads will be concatenated into a single file. Following is a guppy basecaller command which is used to basecall from the rawreads.

```
module load guppy/2.3.1

guppy_basecaller \
        -i ../01_raw_reads \
        -s . \
        --flowcell FLO-MIN106 \
        --kit SQK-LSK109 \
        --qscore_filtering \
        --cpu_threads_per_caller 16
```


This raw fasta file is located in the following directory:

```
02_basecall_pass
├── 5074_test_LSK109_30JAN19-reads-pass.fasta
└── sequencing_summary.txt
```

This will be the fasta file we will be using to demonstrate the assembling of long reads.   

## 3.3 Quality Report
In here we use [nanoplot](https://github.com/wdecoster/NanoPlot) a tool developed to evaluate the statistics of long read data of Oxford Nanopore Technologies and Pacific Biosciences. 

Working directory:
```
long_read_assembly/
├── 03_qc
```
Full slrum script called [nanoplot.sh](long_read_assembly/03_qc/nanoplot.sh) is located in the `03_qc` directory.

```
NanoPlot --summary ../02_basecall_pass/sequencing_summary.txt \
        --loglength \
        -o summary-plots-log-transformed \
        -t 10
```

This will create a summary files and figures which can be found in the output directory specified in the command along with a HTML report. 
![](images/HistogramReadlength.png)  
![](images/LengthvsQualityScatterPlot_dot.png)
![](images/LengthvsQualityScatterPlot_kde.png) 



### 3.3.1 Contaminant screening of long reads

Here we will be using [Centrifuge](https://ccb.jhu.edu/software/centrifuge/) software to classify DNA sequences from microbial samples. It uses a novel indexing scheme which uses a relatively small index to perform a fast classification.

The command used:
```
centrifuge -f \
	-x /isg/shared/databases/centrifuge/b+a+v+h/p_compressed+h+v \
	--report-file report.tsv \
	--quiet \
	--min-hitlen 50 \
	-U ../02_basecall_pass/5074_test_LSK109_30JAN19-reads-pass.fasta
```

Command options:
```
Usage:  centrifuge [options]* -x <cf-idx>  -U <r>  [-S <filename>] [--report-file <report>]

<cf-idx>   Index filename prefix (minus trailing .X.cf)
-U         can be comma-separated file lists
--min-hitlen  minimum length of partial hits (default 22, must be greater than 15)
-p/--threads  number of alignment threads
--quiet       print nothing to stderr except serious errors
```
The complete slurm script is called [centrifuge.sh](long_read_assembly/03_centrifuge/centrifuge.sh) can be found in the 03_centrifuge folder. This will produce a report and classification report based on the contaminants. We will use a python script to remove the contaminates from the fasta file. 

Inital coverage =  112X  

Post centrifuge coverage = 64X

### 3.3.2 Contaminant screening of short reads  

Working directory:  
```
long_read_assembly/
├── ILLUMINA_DATA/
│   ├── 02_Kraken/

```
For this species we also had Illumina short reads. In order to check the contaminate in the reads we used Kraken software. Illumina short reads were initially trimmed using Trimmomatic and before running the contaminant screening using Kraken. 

Kraken command:
```
module load kraken/2.0.8-beta
module load jellyfish/2.2.6

kraken2 -db /isg/shared/databases/kraken2/Standard \
        --paired cat_R1.fq cat_R2.fq \
        --use-names \
        --threads 16 \
        --output kraken_general.out \
        --unclassified-out unclassified#.fastq \
        --classified-out classified#.fastq      \
        --report kraken_report.txt \
        --use-mpa-style
```
Complete slurm script is called [Kraken.sh](long_read_assembly/ILLUMINA_DATA/02_Kraken/Kraken.sh).
This will create the sequence output, report and a summary report. 
```
├── classified_1.fastq
├── classified_2.fastq
├── kraken_report.txt
├── kraken_general.out
├── unlassified_1.fastq
└── unlassified_2.fastq
```

According to the Kraken2 report we can identify the following summay using:
```
grep -v "|" kraken_report.txt
```

which will give :
```
d__Bacteria	197855570
d__Eukaryota	1270391
d__Archaea	62343
d__Viruses	15501
s__unidentified	6
```

Also if you look at the *.err file it will show:
```
254808872 sequences (47531.49 Mbp) processed in 1266.999s (12066.7 Kseq/m, 2250.90 Mbp/m).
  199246490 sequences classified (78.19%)
  55562382 sequences unclassified (21.81%)
```

Initial coverage of Illumina data (post Trimmomatic) = 50X  
After Kraken filtering = 11X  





## 3.4 Assembly  
In this step we will be using two asseblers Flye and Shasta on the base called ONT data.   

Working directory: 
```
long_read_assembly/
├── 04_assembly
```

### 3.4.1 Flye Assembly   
Flye assembler takes data from Pacbio or Oxford Nanopore technologies sequencers and outputs polished contigs. It will repeat graph, that is similar in appearance to the De Bruijn graph. The manner in which this graph is assembled reveals the repeats in the genome allowing for the most accurate assembly.  

Command for running flye assembly:
```
flye --nano-raw ../../03_centrifuge/physcomitrellopsis_africana_rmv_contam.fasta \
        --genome-size 500m \
        --threads 32 \
        --out-dir /UCHC/PublicShare/CBC_Tutorials/Genome_Assembly/long_read_assembly/04_assembly/flye
```

Useage information:  
```
usage: flye 

--nano-raw     ONT raw reads 
--genome-size  estimated genome size 
--out-dir      Output directory
--threads      number of parallel threads
```  

Full slurm script is called [flye.sh](long_read_assembly/04_assembly/flye/flye.sh) can be found in the `04_assembly/flye` directory. This will create he following files inside the flye_assembly folder, and the assembled fasta file is called scaffolds.fasta or assembly.fasta.

```
flye/
├── 00-assembly
├── 10-consensus
├── 20-repeat
├── 21-trestle
├── 30-contigger
├── 40-polishing
├── assembly.fasta
├── assembly.fasta.fai
├── assembly.fasta.mmi
├── assembly_graph.gfa
├── assembly_graph.gv
├── assembly_info.txt
├── params.json
└── scaffolds.fasta
```

### 3.4.1 Flye Assembly with Polish 
In this section we will use flye assembly with 3 polishing iterations.  

```
flye --nano-raw ../../02_basecall_pass/5074_test_LSK109_30JAN19-reads-pass.fasta \
        --genome-size 500m \
        --threads 32 \
        --iterations 3 \
        --out-dir .
```

Flye usage:
```
usage: flye ( --nano-raw ) file1 [file_2 ...] [options]

optional arguments:
--nano-raw path [path ...]   ONT raw reads
--genome-size size           estimated genome size (for example, 5m or 2.6g)
--out-dir path               Output directory
--iterations int             number of polishing iterations
--threads int                number of parallel threads
```

The full slurm script is called [flye.sh](long_read_assembly/04_assembly/flye_polish/flye.sh). 

The output will contain:
```
flye_polish/
├── 00-assembly
├── 10-consensus
├── 20-repeat
├── 21-trestle
├── 30-contigger
├── 40-polishing
├── assembly.fasta
├── assembly_graph.gfa
├── assembly_graph.gv
├── assembly_info.txt
├── flye.log
├── flye.sh
├── params.json
└── scaffolds.fasta --> [full-path]/assembly.fasta
```


### 3.4.2 Shasta Assembly  
[Shasta](https://github.com/chanzuckerberg/shasta) assembler is developed keeping in mind to produce a rapid and accurate assembly and polishing. More information on the shastar assembler can be found over [here](https://www.nature.com/articles/s41587-020-0503-6). Here we are using the Shasta assembler for the ONT reads.  
Working directory is:
```
long_read_assembly/
├── 04_assembly/
│   ├── shasta/
```

Command for running shasta assembly:  
```
shasta --input ../../03_centrifuge/physcomitrellopsis_africana_rmv_contam.fasta \
        --Reads.minReadLength 1000 \
        --memoryMode anonymous \
        --memoryBacking 4K \
        --threads 32
```

Useage information:
```
shasta [options]
--input                Names of input files containing reads
--Reads.minReadLength  Read length cutoff. Reads shorter than this number of bases are discarded  
--memoryMode           Specify whether allocated memory is anonymous or backed by a filesystem   
--memoryBacking        Specify the type of pages used to back memory
--threads              Number of threads

```

Full slurm script is called [shasta.sh](long_read_assembly/04_assembly/shasta/shasta.sh) can be found in the `/long_read_assembly/04_assembly/shasta/` directory. After running the command it will create Assembly.fasta file with other files:  

```
shasta/
├── ShastaRun
│   ├── Assembly.fasta
└── shasta.sh

```

## 3.5 Assembly Evaluation
In this section we will be evaluating the initial assembly created from the above assemblers. 
Working directory:
```
long_read_assembly/
├── 05_initial_assembly_evaluation/
```


### 3.5.1 BUSCO
In here we will evaluate the assemblies using BUSCO. 
Working directory:
```
05_initial_assembly_evaluation/
├── busco/
```
Following commands will be used to evaluate the flye and shasta initial assemblies.
*   flye:
```
busco -i ../04_assembly/flye/assembly.fasta \
        -o busco_flye -l /isg/shared/databases/BUSCO/odb10/viridiplantae_odb10 -m genome
```
Complete slurm script called [busco_flye.sh](long_read_assembly/05_initial_assembly_evaluation/busco/busco_flye.sh) can be found in the busco directory.  

*   shasta:
```
busco -i ../04_assembly/shasta/ShastaRun/Assembly.fasta \
        -o busco_shasta -l /isg/shared/databases/BUSCO/odb10/viridiplantae_odb10 -m genome
```
Complete slurm script called [busco_shasta.sh](long_read_assembly/05_initial_assembly_evaluation/busco/busco_shasta.sh) can be found in the busco directory.  

*   flye with polish:  
```
busco -i ../../04_assembly/flye_polish/assembly.fasta \
        -o busco_flye_polish -l /isg/shared/databases/BUSCO/odb10/viridiplantae_odb10 -m genome
```
Complete slurm script called [busco_flye_polish.sh](long_read_assembly/05_initial_assembly_evaluation/busco/busco_flye_polish.sh) can be found in the busco directory. 


**NOTE:** When running busco you need to copy the augustus config directory to a location which you have the permision to write to and the path to the augustus config directory should be exported beforehand.  

General useage of the command:
```
usage: busco -i [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] -m [MODE] [OTHER OPTIONS] 
```

The command options we will be using:
```
-i FASTA FILE   Input sequence file in FASTA format
-l LINEAGE      Specify the name of the BUSCO lineage
-o OUTPUT       Output folders and files will be labelled with this name
-m MODE         BUSCO analysis mode
					- geno or genome, for genome assemblies (DNA)
					- tran or transcriptome, for transcriptome assemblies (DNA)
					- prot or proteins, for annotated gene sets (protein)
```

Summary of the inital assembly assesment using BUSCO: 
![](images/long_read_inital_assembly_busco_figure.png)



### 3.5.2 QUAST  
[QUAST](http://quast.sourceforge.net/quast) will be used to evaluate genome assemblies. We will be using the program QUAST which will give us the number of contigs, total length and N50 value; the data we are most interested in.

Working directory:
```
05_initial_assembly_evaluation/
├── quast/
```

*  flye assembly :  
```
quast.py ../04_assembly/flye/Assembly.fasta \
        --threads 8 \
        -o quast_flye

```
Complete slurm script called [quast_flye.sh](long_read_assembly/05_initial_assembly_evaluation/quast/quast_flye.sh) can be found in the quast directory.

*  shasta assembly :  
```
quast.py ../../04_assembly/shasta/ShastaRun/Assembly.fasta \
        --threads 8 \
        -o quast_shasta
```
Complete slurm script called [quast_shasta.sh](long_read_assembly/05_initial_assembly_evaluation/quast/quast_shasta.sh) can be found in the quast directory.

*  flye assembly with polish:
```
quast.py ../../04_assembly/flye_polish/assembly.fasta \
        --threads 16 \
        -o quast_flye_polish
```
Complete slurm script called [quast_flye_polish.sh](long_read_assembly/05_initial_assembly_evaluation/quast/quast_flye_polish.sh) can be found in the quast directory.


Summary of the inital assembly assesment using quast:  

|             |  Flye     |  Shasta     |  Flye with Polish     |  
 ------------ |:---------: | :---------: | ---------: 
 contigs (>= 0 bp)    | 9835  | 10134  | 9874  |
 contigs (>= 1000 bp) | 9397  |  7942  | 9562  | 
 contigs (>= 5000 bp) | 7671  |  6278 | 7991   | 
 contigs (>= 10000 bp) | 6395 |  4826 | 6643  |
 contigs (>= 25000 bp) | 4292 |  2463 | 4346  |
 contigs (>= 50000 bp) | 2747  |  904 | 2706   |
 Total length (>= 0 bp)| 539643975 |  220312440 | 507805683 |
 Total length (>= 1000 bp)| 539340131 | 219821137 | 507598179 |
 Total length (>= 5000 bp)| 534538564 | 214920226 | 503143942 |
 Total length (>= 10000 bp)| 525247378 | 204318472 | 493296072 |
 Total length (>= 25000 bp)| 491112867 |  165543884 | 456093461 |
 Total length (>= 50000 bp)| 435500086 |  110696715 | 396803705 |
 **no. of contigs**     |  9806  | 8261 | 9828 |
 Largest contig     | 6909730 |  4103446 | 6326904 |
 Total length       | 539633067 |  220054765 | 507788936 |
 GC (%)      | 38.13  |  43.63 | 38.00 |
 **N50**         | 152392  |  50297 | 128968 |
 N75         | 646489  |  25268 | 56203  |
 L50         |  823 |  891 | 910 |
 L75         |  2207 |  2443 | 2406 |  




## 3.6 Error correction   
Polishing the assemblies using Medaka. By polishing it means calculates a better consensus sequence for a draft genome assembly, find base modifications, and call SNPs with respect to a reference genome.  

Medaka command:   
```
medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR}  -t 16
```

Command options used:
```
medaka_consensus [-h] -i <fastx>

-i  fastx input basecalls (required)
```


*  #### Error correction of flye assembly:
Working directory:
```
long_read_assembly/
├── 06_error_correction
│   ├── flye_assembly/
```

Complete slurm script is called [medaka.sh](long_read_assembly/06_error_correction/flye_assembly/medaka.sh) is located in the `flye_assembly` directory.

This will produce the following files along with the error corrected *consensus.fasta* file:  
```
flye_assembly
├── calls_to_draft.bam
├── calls_to_draft.bam.bai
├── consensus.fasta
├── consensus_probs.hdf
└── medaka.sh
```

*  #### Error correction of shasta assembly: 
``` 
long_read_assembly/
├── 06_error_correction
│   ├── flye_assembly/
```
Complete slurm script is called [medaka.sh](long_read_assembly/06_error_correction/shasta_assembly/medaka.sh) is located in the `shasta_assembly` directory.

This will produce the following files along with the error corrected *consensus.fasta* file:  
```
shasta_assembly/
├── calls_to_draft.bam
├── calls_to_draft.bam.bai
├── consensus.fasta
├── consensus_probs.hdf
└── medaka.sh
```

## 3.7 Evlaluation after error correction   

In this section we will evaluate the assembly after error correction. 

Working directory:  
```
long_read_assembly/
├── 07_evaluation_after_error_correction/
```

### 3.7.1 BUSCO
BUSCO evaluation:  
```
07_evaluation_after_error_correction/
├── busco/
```

*  flye:  
```
export AUGUSTUS_CONFIG_PATH=../../augustus/config

busco -i ../../06_error_correction/flye_assembly/consensus.fasta \
        -o busco_flye -l /isg/shared/databases/BUSCO/odb10/viridiplantae_odb10 -m genome
```
Complete slurm script called [busco_flye.sh](long_read_assembly/07_evaluation_after_error_correction/busco_flye.sh) can be found in the busco directory.  

* shasta:  
```
export AUGUSTUS_CONFIG_PATH=../../augustus/config
busco -i ../../06_error_correction/shasta_assembly/consensus.fasta \
        -o busco_shasta -l /isg/shared/databases/BUSCO/odb10/viridiplantae_odb10 -m genome
```  

Complete slurm script called [busco_shasta.sh](long_read_assembly/07_evaluation_after_error_correction/busco_shasta.sh) cam be found in the busco directory.  

Summary of the statistics:
```
flye assembly stats after error correction:
C:82.4%[S:74.4%,D:8.0%],F:9.4%,M:8.2%,n:425  

shasta assembly stats after error correction:
C:69.4%[S:65.9%,D:3.5%],F:10.8%,M:19.8%,n:425
```

![](images/07_after_error_correction_busco_figure.png)  




## 3.8 Identifying and removing duplicate regions  

Current working directory:  
```
long_read_assembly/
├── 08_purge
```

[Purge Haplotigs](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/) is a pipeline to help with curating genome assemblies. It assures that there is not a combination of sequences between contigs and haplotigs. It uses a system that uses the mapped reads that you assembled and Minimap2 to assess which contigs should be kept in the assembly.

We use this because some parts of a genome may have a very high degree of heterozygosity which causes contigs for both haplotypes of that part of the genome to be assembled as separate primary contigs, rather than as a contig with a haplotig which may cause an error in analysis.

What the purge halpotigs algortihm does is identify pairs of contigs that are syntenic and group them as haplotig. The pipeline uses mapped read coverage and blast/lastz alignments to determine which contigs to keep for the assembly. Dotplots are produced for all flagged contig matches to help the user analyze any remaining ambiguous contigs.

This part must be done in seperate steps as the parameters in each part depend on the results of the previous steps.  

#### Running purge haplotigs: 
When running purge haplotigs you must run the following command separately one after the other. A detailed tutorial on how to run the purge haplotigs can be found [here](https://bitbucket.org/mroachawri/purge_haplotigs/wiki/Tutorial).

##### Pre processing step:
First you need to map your reads to the assembly using minimap2, after that need to sort and index the bam file with samtools index.  

```
minimap2 -t 16 -ax map-ont ${ref} ${read_file}\
       | samtools view -hF 256 - \
	   | samtools sort -@ 16 -m 1G -o ${bam} -T ${tmp.ali} 
```

##### Step 1:
Then using the following command generate a coverage histogram.

```
purge_haplotigs  readhist  -b aligned.bam  -g genome.fasta  [ -t threads ]

REQUIRED:
-b / -bam       BAM file of aligned and sorted reads/subreads to the reference
-g / -genome    Reference FASTA for the BAM file.

OPTIONAL:
-t / -threads   Number of worker threads to use, DEFAULT = 4, MINIMUM = 2
```
This will produce histogram file in png format and a BEDTools ‘genomecov’ output file which we will be needing for the step2:

##### Step 2:
Using the above histogram the cutoff values for low read depth, high read depth cutoff need to be selected. Low read cutoff been selected as the midpoint between the haploid and diploid peaks.
```
purge_haplotigs  contigcov  -i aligned.bam.genecov  -l 30  -m 80  -h 145 

REQUIRED:
-i / -in        The bedtools genomecov output that was produced from 'purge_haplotigs readhist'
-l / -low       The read depth low cutoff (use the histogram to eyeball these cutoffs)
-h / -high      The read depth high cutoff
-m / -mid       The low point between the haploid and diploid peaks

OPTIONAL:
-o / -out       Choose an output file name (CSV format, DEFAULT = coverage_stats.csv)
-j / -junk      Auto-assign contig as "j" (junk) if this percentage or greater of the contig is 
                low/high coverage (DEFAULT = 80, > 100 = don't junk anything)
-s / -suspect   Auto-assign contig as "s" (suspected haplotig) if this percentage or less of the
                contig is diploid level of coverage (DEFAULT = 80)
```

This will produce the following files: 
```
├── coverage_stats.csv
```

##### Step 3:
The following command will run to purge the reads:
```
purge_haplotigs purge -b aligned.bam -g ${ref} -c coverage_stats.csv -d -a 60
```

Command options:
```
REQUIRED:
-g / -genome        Genome assembly in fasta format. Needs to be indexed with samtools faidx.
-c / -coverage      Contig by contig coverage stats csv file from the previous step.
-d / -dotplots      Generate dotplots for manual inspection.
-a / -align_cov     Percent cutoff for identifying a contig as a haplotig. DEFAULT = 70
```

*  ##### flye assembly purge:
working directory:
```
08_purge/
├── flye/
```  
**Step1:**  
```
purge_haplotigs readhist -b flye_aligned.bam -g ${ref} -t 16
``` 

**Step2:** 
![](images/flye_aligned.bam.histogram.png)  
```
purge_haplotigs  contigcov -i flye_aligned.bam.gencov -l 3 -m 55 -h 195
```

**Step3:** 
```
purge_haplotigs purge -b flye_aligned.bam -g ${ref} -c coverage_stats.csv -d -a 60
```

Complete slurm script called [purge_haplotigs.sh](long_read_assembly/08_purge/flye/purge_haplotigs.sh) can be found in the `flye` directory. It will result in a final curated fasta file called *curated.fasta*.


*  ##### shasta assembly purge:   
working directory:
```
08_purge/
├── shasta/
```  

**Step1:**  
```
purge_haplotigs readhist -b shasta_aligned.bam -g ${ref} -t 16
```  

**Step2:**
![](images/shasta_aligned.bam.histogram.png)  
```
purge_haplotigs  contigcov -i shasta_aligned.bam.gencov -l 5 -m 34 -h 195
```

**Step3:**
```
purge_haplotigs purge -b shasta_aligned.bam -g ${ref} -c coverage_stats.csv -d -a 60
```


Complete slurm script called [purge_haplotigs.sh](long_read_assembly/08_purge/shasta/purge_haplotigs.sh) can be found in the `shasta` directory. It will result in a final curated fasta file called *curated.fasta*.  

## 3.9 Final assembly evaluation   

Working directory:
```
long_read_assembly/
├── 09_final_assembly_evaluation/
```

### 3.9.1 BUSCO 
*  ##### BUSCO evaluation for flye assembly
```
busco -i ../../08_purge/flye/curated.fasta \
        -o flye_curated_busco -l /isg/shared/databases/BUSCO/odb10/viridiplantae_odb10 -m genome
```
complete slurm script called [busco_flye.sh](long_read_assembly/09_final_assembly_evaluation/busco/busco_flye.sh) can be found in the `busco` directory.

The summary of the out will contain in:
```
flye_curated_busco/
└── short_summary.specific.viridiplantae_odb10.flye_curated_busco.txt
```

This will contain the short summary of the evaluation which contains:
```
C:80.2%[S:72.0%,D:8.2%],F:10.6%,M:9.2%,n:425
341     Complete BUSCOs (C)
306     Complete and single-copy BUSCOs (S)
35      Complete and duplicated BUSCOs (D)
45      Fragmented BUSCOs (F)
39      Missing BUSCOs (M)
425     Total BUSCO groups searched
```



*  ##### BUSCO evaluation for shasta assembly  
```
busco -i ../../08_purge/shasta/curated.fasta \
        -o shasta_curated_busco -l /isg/shared/databases/BUSCO/odb10/viridiplantae_odb10 -m genome
```
complete slurm script called [busco_shasta.sh](long_read_assembly/09_final_assembly_evaluation/busco/busco_shasta.sh) can be found in the `busco` directory.

The summary of the out will contain in:
```
shasta_curated_busco/
└── short_summary.specific.viridiplantae_odb10.shasta_curated_busco.txt
```
This will contain the short summary of the evaluation which contains:
```
C:68.5%[S:65.2%,D:3.3%],F:11.5%,M:20.0%,n:425
291     Complete BUSCOs (C)
277     Complete and single-copy BUSCOs (S)
14      Complete and duplicated BUSCOs (D)
49      Fragmented BUSCOs (F)
85      Missing BUSCOs (M)
425     Total BUSCO groups searched
```  

![](images/09_final_assembly_evaluation_busco_figure.png)


### 3.9.2 QUAST  
QUAST will be used to evaluate genome assemblies

*  ##### QUAST evaluation for flye assembly  

```
quast.py ../../08_purge/flye/curated.fasta \
        --threads 8 \
        -o quast_flye
```
Complete slurm script called [quast_flye.sh](long_read_assembly/09_final_assembly_evaluation/quast/quast_flye.sh) can be found in the quast directory.

Summary of the statistics can be found in:
```
quast_flye
├── report.txt
```

*  ##### QUAST evaluation for shasta assembly 

```
quast.py ../../08_purge/shasta/curated.fasta \
        --threads 8 \
        -o quast_shasta
```
Complete slurm script called [quast_shasta.sh](long_read_assembly/09_final_assembly_evaluation/quast/quast_shasta.sh) can be found in the quast directory.

Summary of the statistics can be found in:
```
quast_shasta/
├── report.txt
```

Statistics summary of flye and shasta assemblies:

|             |  flye     |  shasta     |      
 ------------ |:---------: | :---------: |
 contigs (>= 0 bp)    | 6549  | 8119  |
 contigs (>= 1000 bp) | 6392   |  7344 | 
 contigs (>= 5000 bp) | 5919   |  6005 |  
 contigs (>= 10000 bp) | 5366  |  4727 |
 contigs (>= 25000 bp) | 3971   |  2437  |
 contigs (>= 50000 bp) | 2611    |  887 |
 Total length (>= 0 bp)| 492754620 |  212410761 | 
 Total length (>= 1000 bp)| 492661682 | 212090599 | 
 Total length (>= 5000 bp)| 491351568 | 208075243|
 Total length (>= 10000 bp)| 487138762| 198692220 |
 Total length (>= 25000 bp)| 463778615 |  161021099 |
 Total length (>= 50000 bp)| 414473760 |  106478172 |
 **no. of contigs**     |  6501  | 7587 | 
 Largest contig     | 6892258 |  4089676 |
 Total length       | 492739170 |  212267457 |
 GC (%)      | 38.43  |  43.52 | 
 **N50**         | 160624  |  50177 | 
 N75         | 72963  |  25781 |
 L50         |  719 |  881|
 L75         |  1872|  2366|


### 3.9.3 Minimap2  

In here we will be aligning our reads to the assembed genome to see the alignment rate. Mapping of the reads will be done using the minimap2 and the generated SAM file then will be sorted into BAM format file. 
```
minimap2 -t 16 -ax map-ont ${ref} ${read_file} \
        | samtools sort -@ 16 -m 2G -o ${aligned.bam} -T ${tmp.ali}

bamtools stats -in  ${aligned.bam}

```

*  Complete slrum script for aligning the reads to the flye assembly is called [minimap2_flye.sh](long_read_assembly/09_final_assembly_evaluation/minimap2/minimap2_flye.sh) can be found in the minimap2 directory.  

*  Complete slrum script for aligning the reads to the flye assembly is called [minimap2_shasta.sh](long_read_assembly/09_final_assembly_evaluation/minimap2/minimap2_shasta.sh) can be found in the minimap2 directory. 

|             |  flye     |  shasta     |      
 ------------ |:---------: | :---------: |
 Total reads    | 25964212  | 26949492  |
 Mapped reads | 22520313 (86.736%)   |  18231482 (67.6506%) | 
 Forward strand | 14787262 (56.9525%)   |  17814705 (66.104%) |  
 Reverse strand | 11176950 (43.0475%)  |  9134787 (33.896%) |


## 4. Hybrid Genome Assembly

### Introduction
In this section we will be working with hybrid assemblers which will be compatible with long read PacBio Data and Short read Illumina data. Nanopore and PacBio are currently both the main long read sequencing technologies but the major differences in them are that PacBio reads a molecule multiple times to generate high-quality consensus data while Nanopore can only sequence a molecule twice. As a result, PacBio generates data with lower error rates compared to Oxford Nanopore.

To perform a hybrid assembly it requires you have both short and long read data to complete the genome. Hybrid assembly uses short read data to resolve ambiguity in the long read data as it is assembled. For this tutorial we are using data from a boxelder (Acer_negundo) genome. 

##### Pacbio Data
First step is the conversion of RS-II pacbio data into fasta format, which is done using [pbh5tools](https://github.com/PacificBiosciences/pbh5tools/blob/master/doc/index.rst) provided by pacificbiosciences. `bash5tools.py` can extract read sequences and quality score values from raw and circular consensus sequencing reads to create fastq and fasta files.
```
bash5tools.py --outFilePrefix [PREFIX] --outType fasta --readType [Type] input.bas.h5  

input.bas.h5          input .bas.h5 filename
--outFilePrefix OUTFILEPREFIX
--readType {ccs,subreads,unrolled}
--outType OUTTYPE     output file type (fasta, fastq) [fasta]  
```

Our Pacbio data is located in the `Pacbio_DATA` directory
```
Pacbio_DATA/
└── acne_pb.fasta
```

##### Preprocessing with PacBio CCS  
CCS takes multiple subreads of the same molecule and combines them using a statistical model to produce one accurate consensus sequence (HiFi read), with base quality values. For more information, refer to the [PacBio](https://www.pacb.com/smrt-science/smrt-sequencing/hifi-reads-for-highly-accurate-long-read-sequencing/) Website. For more information on CCS please refer to [this](https://www.pacb.com/videos/tutorial-circular-consensus-sequence-analysis-application/) PacBio CCS Tutorial. And if you have a UConn PacBio account please refer to [this](https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/pacbio-v7/) tutorial.  



##### Illumina Data
The short reads that we are going to use in this assembly is located in the `ILLUMINA_DATA/` directory:
```
ILLUMINA_DATA/
├── XUMNS_20180703_K00134_IL100105003_S37_L007_R1.fastq.gz
├── XUMNS_20180703_K00134_IL100105003_S37_L007_R2.fastq.gz
├── XUMNS_20180703_K00134_IL100105003_S37_L008_R1.fastq.gz
└── XUMNS_20180703_K00134_IL100105003_S37_L008_R2.fastq.gz
```

## Assembly

## A. Short Read Assembly  

### A.1 Illumina Short Read Assembly using Masurca
In here we will assemble the short read illumina reads using masurca.
Working directory:
```
hybrid_assembly/
├── Acer_negundo/
│   └── Short_Read_Assembly/
│       └── 01_masurca_assembly/
``` 

config file will look like:  
```
DATA
PE = aa 320 20 ../../ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L007_R1.fastq.gz ../../ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L007_R2.fastq.gz
PE = ab 320 20 ../../ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L008_R1.fastq.gz ../../ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L008_R2.fastq.gz
END

PARAMETERS
EXTEND_JUMP_READS=0
GRAPH_KMER_SIZE = auto
USE_LINKING_MATES = 0
USE_GRID=0
GRID_ENGINE=SLURM
GRID_QUEUE=general
GRID_BATCH_SIZE=100000000
LHE_COVERAGE=25
MEGA_READS_ONE_PASS=0
LIMIT_JUMP_COVERAGE = 300
CA_PARAMETERS =  cgwErrorRate=0.25
CLOSE_GAPS=1
NUM_THREADS = 32
JF_SIZE = 4500000000
SOAP_ASSEMBLY=0
FLYE_ASSEMBLY=0
END

```

The command we use is:
```
masurca config.txt 

./assemble.sh
```
The slurm script called [masurca_assembly.sh](hybrid_assembly/Acer_negundo/Short_Read_Assembly/01_masurca_assembly/masurca_assembly.sh) can be found in *01_masurca_assembly* directory.  

This will create the assembly using only the Illumina reads.  
Final assembly file will be in the **CA** directory:
```
01_masurca_assembly/
├── CA/
│   ├── final.genome.scf.fasta
```

###  A.2 Assembly Evaluation   

Working directory:  
```
Short_Read_Assembly/
├── 02_evaluation/
```

#### Quast evaluation  

Final assembly evaluation was done using QUAST:
```
quast.py ../01_masurca_assembly/CA/final.genome.scf.fasta \
        --threads 16 \
        -o masurca_assembly
```

The full slurm script called [quast.sh](hybrid_assembly/Acer_negundo/Short_Read_Assembly/02_evaluation/quast.sh), can be found in the *02_evaluation* directory.   

Resulting statistics on the assembly will be written to the *report.txt* file.  
```
masurca_assembly/
├── report.txt
```

|             |  masurca     |
------------ |:---------: |  
Total length (>= 50000 bp)  |  51091  |  
contigs  |  162899   |  
Largest contig  |  51091  |  
Total length  |  426302651  |  
N50   |  4563  |  

  
#### BUSCO evaluation  
BUSCO evaluation:  
```
busco -i  ../01_masurca_assembly/CA/final.genome.scf.fasta \
        -o busco -l /isg/shared/databases/BUSCO/odb10/viridiplantae_odb10 -m genome
```

The full slurm script is called [busco.sh](hybrid_assembly/Acer_negundo/Short_Read_Assembly/02_evaluation/busco.sh). 


The summary of the out will contain in:
```
busco/
└── short_summary.specific.viridiplantae_odb10.busco.txt
```
This will contain the short summary of the evaluation which contains:
```
C:79.0%[S:77.4%,D:1.6%],F:18.4%,M:2.6%,n:425
336     Complete BUSCOs (C)
329     Complete and single-copy BUSCOs (S)
7       Complete and duplicated BUSCOs (D)
78      Fragmented BUSCOs (F)
11      Missing BUSCOs (M)
425     Total BUSCO groups searched
```  



## B. Long Read Assembly  
Working directory:

```
hybrid_assembly/
├── Acer_negundo/
│   └── Long_Read_Assembly/
```

### B.1 Pacbio long read assembly using flye  
In this section we will assemble the pacbio long reads using flye assembler.  
Working directory:
```
Long_Read_Assembly/
├── 01_flye_assembly/
```

Command used:  
```
flye --pacbio-raw ../Pacbio_DATA/acne_pb.fasta \
        --genome-size 300m \
        --threads 32 \
        --out-dir . 
```

The full slurm script called [flye.sh](hybrid_assembly/Acer_negundo/Long_Read_Assembly/01_flye_assembly/flye.sh) cam be found in *01_flye_assembly* directory.  


#### Pacbio Assembly using Canu
In here we will be only using pacbio reads to do the assembly. 

```
module load gnuplot/5.2.2
module load canu/2.1.1

canu useGrid=true \
        -p Acer_negundo -d canu_out \
        genomeSize=4.5g \
        -pacbio /UCHC/PublicShare/CBC_Tutorials/Genome_Assembly/hybrid_assembly/Acer_negundo/Pacbio_DATA/acne_pb.fasta gridOptions="--partition=himem --qos=himem  --mem-per-cpu=8000m --cpus-per-task=24"
```  


### B.2 Evaluation  
Now we will evaluate the final assembly created by flye assembler.  

Working directory:  
```
Long_Read_Assembly/
├── 02_evaluation/
```  

#### Quast evaluation  
Quast command used: 
```
quast.py ../01_flye_assembly/assembly.fasta \
        --threads 16 \
        -o quast 
```

The full slurm script called [quast.sh](hybrid_assembly/Acer_negundo/Long_Read_Assembly/02_evaluation/quast.sh) can be found in *02_evaluation* directory.  

Resulting statistics on the assembly will be written to the *report.txt* file.  
```
masurca_assembly/
├── report.txt
```

|             |  flye     |
------------ |:---------: |  
Total length (>= 50000 bp)  |  441362208  |  
contigs  |  1684  |  
Largest contig  |  8254085  |  
Total length  |  4451717887  |  
N50   |  1762694  |  

#### BUSCO evaluation  
BUSCO evaluation:  
```
busco -i  ../01_flye_assembly/assembly.fasta \
        -o busco -l /isg/shared/databases/BUSCO/odb10/viridiplantae_odb10 -m genome
```

The full slurm script is called [busco.sh](hybrid_assembly/Acer_negundo/Long_Read_Assembly/02_evaluation/busco.sh). 


The summary of the out will contain in:
```
busco/
└── short_summary.specific.viridiplantae_odb10.busco.txt
```
This will contain the short summary of the evaluation which contains:
```
C:98.1%[S:95.5%,D:2.6%],F:1.2%,M:0.7%,n:425
417     Complete BUSCOs (C)
406     Complete and single-copy BUSCOs (S)
11       Complete and duplicated BUSCOs (D)
5      Fragmented BUSCOs (F)
3       Missing BUSCOs (M)
425     Total BUSCO groups searched
```


## C. Hybrid assembly  

Working directory:

```
hybrid_assembly/
├── Acer_negundo/
│   └── Hybrid_Assembly/
```

### C.1  Assembly using Masurca  

Working directory:
```
Hybrid_Assembly/
├── 01_masurca_assembly/
```

We will be using Illumina reads and long read pacbio reads to construct a hybrid assembly using masurca. For masurca assembly we need a [configuration file](hybrid_assembly/masurca_assembly/config.txt) which directs to the short read and long reads. 

```
DATA
PE = aa 320 20 ../../ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L007_R1.fastq.gz ../../ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L007_R2.fastq.gz
PE = ab 320 20 ../../ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L008_R1.fastq.gz ../../ILLUMINA_DATA/XUMNS_20180703_K00134_IL100105003_S37_L008_R2.fastq.gz
PACBIO = /UCHC/PublicShare/CBC_Tutorials/Genome_Assembly/hybrid_assembly/Acer_negundo/Pacbio_DATA/acne_pb.fasta 
END

PARAMETERS
EXTEND_JUMP_READS=0
GRAPH_KMER_SIZE = auto
USE_LINKING_MATES = 0
USE_GRID=0
GRID_ENGINE=SLURM
GRID_QUEUE=himem
GRID_BATCH_SIZE=100000000
LHE_COVERAGE=25
MEGA_READS_ONE_PASS=0
LIMIT_JUMP_COVERAGE = 300
CA_PARAMETERS =  cgwErrorRate=0.25
CLOSE_GAPS=1
NUM_THREADS = 32
JF_SIZE = 4500000000
SOAP_ASSEMBLY=0
FLYE_ASSEMBLY=0
END
```

The full slurm script called [masurca.sh](hybrid_assembly/masurca_assembly/masurca.sh) is located in the `masurca_assembly` directory.

Final assembly file will be in the CA directory:
```
01_masurca_assembly/
├── CA.mr.41.17.20.0.02/
│   ├── final.genome.scf.fasta
```

### C.2 Evaluation   

Working directory:
```
Hybrid_Assembly/
├── 02_evaluation/
```

#### Quast evaluation  
Quast command used: 
```
quast.py ../01_masurca_assembly/CA.mr.41.17.20.0.02/final.genome.scf.fasta \
        --threads 16 \
        -o masurca_assembly
```
  

he full slurm script called [quast.sh](hybrid_assembly/Acer_negundo/Hybrid_Assembly/02_evaluation/quast.sh) can be found in *02_evaluation* directory.  

Resulting statistics on the assembly will be written to the *report.txt* file.  
```
masurca_assembly/
├── report.txt
```  

|             |  masurca     |
------------ |:---------: |  
Total length (>= 50000 bp)  |  485855121  |  
contigs  |  1993  |  
Largest contig  |  6444170  |  
Total length  |  509195518  |  
N50   |  143  |  


#### BUSCO evaluation  
BUSCO evaluation:  
```
busco -i  ../01_masurca_assembly/CA.mr.41.17.20.0.02/final.genome.scf.fasta \
        -o busco -l /isg/shared/databases/BUSCO/odb10/viridiplantae_odb10 -m genome
```

The full slurm script is called [busco.sh](hybrid_assembly/Acer_negundo/Hybrid_Assembly/02_evaluation/busco.sh). 


The summary of the out will contain in:
```
busco/
└── short_summary.specific.viridiplantae_odb10.busco.txt
```
This will contain the short summary of the evaluation which contains:
```
C:98.6%[S:90.6%,D:8.0%],F:0.7%,M:0.7%,n:425
419     Complete BUSCOs (C)
385     Complete and single-copy BUSCOs (S)
34       Complete and duplicated BUSCOs (D)
3      Fragmented BUSCOs (F)
3       Missing BUSCOs (M)
425     Total BUSCO groups searched
```



### References
*   NanoPack: visualizing and processing long-read sequencing data, Bioinformatics, Volume 34, Issue 15, 01 August 2018, Pages 2666–2669. 
*   Shasta: Nanopore sequencing and the Shasta toolkit enable efficient de novo assembly of eleven human genomes, Nature Biotechnology volume 38, pages1044–1053(2020)