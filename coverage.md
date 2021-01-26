#### How to calculate the coverage   

```
Coverage    =  Total number of bases / genome size

            =  (read count * read length) / genome size
```

So let's see how can we calculate the total number of bases in our data.  

What we have here is paired end data in the form of a fastq file. Each sequence of the FASTQ file record is associated with 4 lines:  
```
@M00704:48:000000000-AH6H7:1:1101:15536:1333 1:N:0:2 
ACATATACTTCTTCATCAGATTCATTATTATCAGTGAATGATGCTTGATCTTGTTCAGACCCCT
+
1>>1>FFFDDFDG3FGB3F1EGH3FEGFFDBFFFHFHHHFGDGHFD1BGHFHH1FGH2GHHHGG
```   

There are multiple ways to calculate the number of sequence in a FASTQ file and in here I am showing one of such methods to calcuate the number of sequences in the FASTQ file. So to calculate the number of sequences in a FASTQ file you need to calculate the total number of lines and then divide it by 4.   

In here I am using a awk command to calculate the number of sequences in a FASTQ file:   
```
awk '{s++}END{print s/4}' Sample_R1.fastq
```  

The above command will give you the number of sequences in a single file. As we are using paired end data, we need to multipy that number by 2 to get the total number of reads in the library. 

```
Total number of reads        = (number of reads in a single file) * 2   
```

As the readlength of a read is 250bp, to find the total number of bases, we need to multipy the Total number of reads by read length of 250bp.

```
Total number of bases       = (Total number of reads) * 250  
```  

As now we now know the total number of bases, and the estimated genome sieze, we can now calculate the **coverage** using the above values.

