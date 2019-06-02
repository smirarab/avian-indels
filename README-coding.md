Scripts to perform indel coding

This code should run in almost all *nix environments. We ran it on a Mac OS X (several
versions, including v. 10.13). The goal of this code is to code indels and convert the 
data into a format that is easy to use for analyses.

The actual indel coding was performed by 2matrix, a perl program available from 
https://github.com/nrsalinas/2matrix

The original code used 2xread, a precursor to 2matrix. This is why some of the pipeline
refers to 2xread. It also assumes the directory structure from the Jarvis et al. data,
which can be downloaded from:

ftp://parrot.genomics.cn/gigadb/pub/10.5524/101001_102000/101041/FASTA_files_of_loci_datasets/

NOTE that the Jarvis data is described in Jarvis, E. D., Mirarab, S., Aberer, A. J., Li, B., 
Houde, P., Li, C., et al. (2015). Phylogenomic analyses data of the avian phylogenomics 
project. GigaScience, 4(1), 4.

### Programs:

#### introns:

* code-individual-introns.pl
* extract-intron.cpp
* jarvis-taxa-authority.txt
* multi-indel-code.sh
* reformat-2xread-output.sh

`extract-intron.cpp` is straightforward to compile. Simply use:

```
g++ -c extract-intron.cpp
g++ -o extract-intron extract-intron.o
```

(I have assumed g++ but the program should be robust)

NOTE that this program assumes a partition file defining individual introns is available
(this is true for the Jarvis intronic data). It will break the data up into individual
introns, indel code those introns, and then produce concatenated files that include all
indels in a locus. Those files include both relaxed phylip format and nexus format files;
the nexus files include partitions for the individual introns and for indels of different
sizes.

A typical sets block for a nexus file looks like this:

```
 Begin sets;
 	charset locus_5_intron12 = 1-379 ;
 	charset locus_5_intron13 = 380-997 ;
 	charset locus_5_intron14 = 998-1335 ;
 	charset locus_5_intron15 = 1336-1665 ;
 	charset locus_5_intron16 = 1666-2103 ;
 	[ four charsets based on indel length: ]
 	[   charset len2plus   = >1 bp indels  ]
 	[   charset len10plus  = >9 bp indels  ]
 	[   charset len50plus  = >49 bp indels ]
 	[   charset len100plus = >99 bp indels ]
 	charset locus_5_len2plus = 1 3 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21 22 ... etc ... ;
 	charset locus_5_len10plus = 10 13 14 15 16 17 18 19 20 21 22 29 33 43 49 53 55 ... etc ...  ;
 	charset locus_5_len50plus = 14 108 161 296 312 313 361 445 467 507 569 632 639 ... etc ... ;
 	charset locus_5_len100plus = 14 108 161 312 467 507 639 661 676 687 844 845 909 ... etc ...  ;
 End;
```

Log files with information about the indels are also generated.

To code all introns run multi-indel-code.smulti-indel-code.sh like this:

```
./multi-indel-code.sh 1000 10002 10006 10007 10008 10009 1001 10011 ... etc
```

(the numbers of are the intron locus numbers from Jarvis et al.)


####  UCEs:

* code-individual-uce.pl
* jarvis-taxa-authority.txt
* reformat-2xread-output.sh
* uce_indel_code.sh

Unlike introns UCEs correspond to a contiguous region. The shell script to code all introns
(uce_indel_code.sh) is hard wired to run code-individual-uce.pl on each locus

