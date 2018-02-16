This is a workflow for cleaning Illumina paired end 150 metagenomic shotgun reads and 
assembling them into contigs. The contigs were then used in a Hi-C cluster analysis but 
this cleaning and assembly pipeline can be used whenever cleaned or assembled reads are 
needed for any analysis.

Project website and data repository: [Using Hi-C to Track Plasmids and Antibiotic
Resistance in Microbial Communities](https://osf.io/gr2d7/)

By Elle J. Kohler, kohl5779@vandals.uidaho.edu

Project paper: link and citation coming soon

## Read Cleaning
The sequenced metagenomic shotgun reads were processed/cleaned prior to assembling them 
into contigs. The cleaning process consisted of three steps: removing duplicates introduced
 by PCR amplification, merging overlapping reads, and trimming off low quality read ends. 
 All programs were run via command line on a remote server in order to have the computational 
 power to run them quickly. The commands used to run each program are given as 
 well as some sample outputs for [Community 4a](https://osf.io/zacf7/). The metagenomic shotgun 
 reads for Community 4a can be downloaded [here](https://osf.io/pvj64/) and used to run 
 this tutorial and duplicate the results. Remember that depending on what directory you 
 are running your analysis from, paths to the files may need to be specified in your commands.
 
### 1. Remove Duplicates
The metagenomic shotgun reads were first sorted using Super-Deduper to remove PCR duplicates 
introduced by the PCR amplification done prior to sequencing. 
[Super-Deduper](https://github.com/dstreett/Super-Deduper) is an open source application 
and was used with default parameters on the raw, paired-end reads via the command:

```> super_deduper -1 AD004_S2_R1_001.fastq.gz -2 AD004_S2_R2_001.fastq.gz```

This will produce the files: **output_nodup_PE1.fastq** and **output_nodup_PE2.fastq**, 
as well as the following output to screen:

```Final:| reads: 25496938 | duplicates: 1971431 | reads_written: 23525507 | percent: 7.73 | discarded: 8593 | total_seconds: 334.63 | reads/sec: 76194.42```

This shows how many read pairs were present, 25496938, and how many duplicates were 
discarded, 1971431. If the discarded reads were a larger percentage of the total reads it 
could suggest that the sample had been over sequenced for its level of diversity.

### 2. Merge Overlapping Reads
The output from Super-Deduper was used as input for [Flash2](https://github.com/dstreett/FLASH2). 
Flash2 merges paired-end reads that were sequenced from DNA fragments shorter than the 
combined length of the reads. For example, when a 250bp sequence fragment is sequenced 
using a paired end 150 kit (as was used for the sample dataset), each end is sequenced for 
150bp. The resulting reads will then overlap by 25bp. Merging these reads prior to metagenomic 
assembly results in a simpler dataset to help the assembly take less time and computational power. 
Flash2 was used with default parameters plus several additional options: -M 200 to set the 
maximum overlap length to 200 (overkill since our reads were only expected to be ~150bp in 
each direction, can be adjusted if longer or shorter reads are being used), -O checks for 
read overlaps at both ends of the reads, -Q 20 sets the quality score (reads are cut off 
if they fall below this value, but quality score is only used to break 
ties in the case of multiple overlaps), -C 70 sets the percentage of the read to be cut off if 
it falls below the quality score.

```> flash2 -M 200 -O -C 70 -Q 20 output_nodup_PE1.fastq output_nodup_PE2.fastq```

These were slightly aggressive cleaning parameters but they were chosen since many of the samples 
were sequenced quite deeply for their complexity and we could afford to throw out any low 
quality reads. Flash2 output three files: 
**out.notCombined_1.fastq**, **out.notCombined_2.fastq**, and **out.extendedFrags.fastq**. 
The file **out.extendedFrags.fastq** contained the reads that were overlapped and thus no 
longer had a pair. Flash2 also produced the following output to screen at the end of the run:

```[FLASH] Read combination statistics:
[FLASH]     Total pairs:       23525507
[FLASH]     Discarded pairs:   2705
[FLASH]     Percent Discarded: 0.01%
[FLASH]     Combined pairs:    5553696
[FLASH]         Innie pairs:   3917323 (70.54% of combined)
[FLASH]         Outie pairs:   1636373 (29.46% of combined)
[FLASH]     Uncombined pairs:  17969106
[FLASH]     Percent combined:  23.61%
```
The most important statistic here is the percentage of reads that overlapped and were 
combined: 23.61%. If this number is very high it could suggest that the sequence run ran 
longer than necessary for the size of the fragments being sequenced. Flash2 also discarded 
a very small percentage (0.01%) due to low quality in overlapping regions.

### 3. Trim Low Quality Ends
[Sickle](https://github.com/najoshi/sickle) was used to trim both the overlapped and 
non-overlapped reads output by flash2. Reads produced from most sequencing technologies 
have progressively lower quality approaching the 3'-ends of the reads and sometimes the 5'-ends as well. 
Trimming off these low quality bases improves the read quality and aids assembly later. 
Sickle will also discard reads based upon a length threshold if specified. This means that 
reads that have been quality trimmed to too short of a length can be discarded. This 
capability was utilized because if a read is too short it only makes for a larger dataset 
without adding much useful genetic information and could slow down assembly. Sickle was 
used with the following parameters: -n removes all reads with an N in them since this 
denotes low quality (this again is a somewhat aggressive parameter that we could afford to
use with our data), --length-threshold 75 discards reads that were trimmed to a length
shorter than 75bp, --qual-threshold 20 sets the quality threshold for trimming, --qual-type 
sanger is used for reads processed using CASAVA 1.8 or higher as was the case with our 
modern, Illumina reads. The paired end reads and single overlapped fragments output by flash2
got processed separately with Sickle. First the paired end reads were processed using the sickle pe command:

```> sickle pe -n --length-threshold 75 --qual-threshold 20 --qual-type sanger -f out.notCombined_1.fastq -r out.notCombined_2.fastq -o cleaned_PE1.fastq -p cleaned_PE2.fastq -s cleaned_SE1.fastq```

The output is three files: **cleaned_PE1.fastq**, **cleaned_PE2.fastq**, and **cleaned_SE1.fastq**. 
There are two files for the forward and reverse reads and one for single reads. 
The single reads were created when the quality of one read in the pair was too low causing 
it to be discarded. The screen output is:

```
PE forward file: out.notCombined_1.fastq 
PE reverse file: out.notCombined_2.fastq

Total input FastQ records: 35938212 (17969106 pairs)

FastQ paired records kept: 33546492 (16773246 pairs)
FastQ single records kept: 1044691 (from PE1: 855432, from PE2: 189259)
FastQ paired records discarded: 302338 (151169 pairs)
FastQ single records discarded: 1044691 (from PE1: 189259, from PE2: 855432)

PE1 Base pairs left removed: 1351646
PE1 Base pairs right removed: 97791786
PE2 Base pairs left removed: 2069799
PE2 Base pairs right removed: 109162395
```
These statistics show how many reads were completely discarded and how many base pairs were
removed from each end of the input reads. They show how many more base pairs were removed 
from the right ends of the reads than the left.

Sickle se was run on the file of single overlapped reads output by flash2. It was run with
the same parameters as sickle pe:

```> sickle se -n --length-threshold 75 --qual-threshold 20 --qual-type sanger --fastq-file out.extendedFrags.fastq --output-file cleaned_SE2.fastq```
 
It output the file: **cleaned_SE2.fastq**, as well as the screen output:

```
SE input file: out.extendedFrags.fastq

Total FastQ records: 5553696
FastQ records kept: 5250686
FastQ records discarded: 303010

Base pairs left removed: 155524
Base pairs right removed: 7594357
```
Again more base pairs were quality trimmed off the right end of the reads than off the left.
The output files, **cleaned_SE1.fastq** and **cleaned_SE2.fastq** were concatenated into 
one file of unpaired reads in preparation for assembly:

```> cat cleaned_SE1.fastq cleaned_SE2.fastq > cleaned_SE.fastq```

## Metagenomic Assembly
The cleaned reads were assembled into longer contiguous sequences (contigs) using the 
assembler [SPAdes version 3.9.0](http://cab.spbu.ru/software/spades/). It was run in meta
mode with the default k-mer sizes of 21, 33, and 55.

```> spades.py --meta -1 cleaned_PE1.fastq -2 cleaned_PE2.fastq -s cleaned_SE.fastq -o spades_output/```

The default parameters allowed it to use up to 40 threads and 400Gb of memory and it took 
~5.5 hours to run. More complex datasets can require a higher amount of allowable memory 
which must be set manually. The output gets put in a directory called **spades_output**. 
This includes a variety of files, the 
two most important of which are: **contigs.fasta** and **spades.log**. The file 
**spades.log** gives information about the run including how long the run took and whether
the run failed or finished. The commands:

```
> head -n 2 contigs.fasta
 
>NODE_1_length_648482_cov_282.424
CCGAAGCCCGCCGATGCGGGCTTCGGCCGTTCCGGGCCCGCCTTTTCGGCGGGCCTTTGC
```
and:

```
> tail -n 2 contigs.fasta 

>NODE_3179_length_56_cov_607
CGGGGTGTGTAGGGCGAATAACGCCATGGGCGTTATCCGCCGATGTCGCGGATGAG
```

will give the outputs shown and can be used to get a quick idea of how well the assembly 
worked. The reads from Community 4a
assembled into 3179 total contigs, the longest of which was 648482bp. The shortest contigs 
were only 56bp and not useful for further analysis. Only those contigs longer than 500bp
were kept for further analysis.

### Assembly Statistics
The contigs longer than 500bp could be separated into a separate file in multiple ways. 
I used the executable script contiglength.py (note that this only works for contig files in
the output format used by SPAdes) to create a list containing the length of the each contig 
using the command:

```> python2.7 contiglength.py contigs.fasta > length4a.txt```

I used the file **length4a.txt** to compute some basic statistics on the assembly in R. 
The script contigstats.R can be used to do this. It will calculate the median and mean 
contig lengths that were assembled (the higher these numbers are, the better the assembly) 
as well as the total number of base pairs assembled (want this number to be close to the 
combined length of each of the replicons present in the sample. For Community 4a, this 
number would be ~17.29Mb and we recovered 17.296Mb. This expected value would not be known 
for uncharacterized samples however.). 
Scrolling through the list of contig lengths in R also makes it easy to find how many contigs 
are longer than 500bp. For Community 4a, 474 of the contigs were longer than 500bp. These 
contigs can be saved into a separate file using the command:

```> sed '/>NODE_475_/Q' contigs.fasta > longcontigs.fasta```

This command is again unique to the contig output format used by SPAdes but could be easily 
adjusted to other output formats. It takes all of the contigs before contig (NODE) 475 and 
copies them into a new file called **longcontigs.fasta**. This file can be used for further 
downstream analysis. If running a Hi-C analysis, save this file and use it for the next 
step in the pipeline which is [cleaning the Hi-C reads](https://osf.io/7n8rx/wiki/home/) 
in preparation for contig clustering based on Hi-C linkages.
