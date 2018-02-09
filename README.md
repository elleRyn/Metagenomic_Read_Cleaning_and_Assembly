### Read Cleaning
The sequenced metagenomic shotgun reads were processed/cleaned prior to assembling them 
into contigs. The cleaning process consisted of three steps: removing duplicates introduced
 by PCR amplification, merging overlapping reads, and trimming off low quality read ends. 
 All programs were run via command line. The commands used to run each tool are given as 
 well as sample outputs for [Community 4a](https://osf.io/zacf7/). The metagenomic shotgun 
 reads for Community 4a can be downloaded [here](https://osf.io/pvj64/) and used to 
 duplicate the results.
#### Remove Duplicates
The metagenomic shotgun reads were sorted using Super-Deduper to remove PCR duplicates 
introduced by the PCR amplification required prior to sequencing. 
[Super-Deduper](https://github.com/dstreett/Super-Deduper) is an open source application 
and was used with default parameters on the raw, paired-end reads via the command:

```> super_deduper -1 AD004_S2_R1_001.fastq.gz -2 AD004_S2_R2_001.fastq.gz```

This will produce the files: **output_nodup_PE1.fastq** and **output_nodup_PE2.fastq**, 
as well as the following output to screen:

```Final:| reads: 25496938 | duplicates: 1971431 | reads_written: 23525507 | percent: 7.73 
| discarded: 8593 | total_seconds: 334.63 | reads/sec: 76194.42```

This shows how many read pairs were present, 25496938, and how many duplicates were 
discarded, 1971431.
#### Merge Overlapping Reads
The output from Super-Deduper was used as input for [Flash2](https://github.com/dstreett/FLASH2). 
Flash2 merges paired-end reads that were sequenced from DNA fragments shorter than the 
combined length of the reads. For example, when a 250bp sequence fragment is sequenced 
using a paired end 150 kit, each end is sequenced for 150bp. The resulting reads then 
overlap by 25bp. Merging these reads prior to metagenomic assembly results in a simpler 
dataset for the assembler to help the assembly take less time and computational power. 
Flash2 was used with default parameters plus several additional options: -M 200 to set the 
maximum overlap length to 200 (overkill since our reads were only expected to be ~150bp in 
each direction), -O checks for read overlaps at both ends of the reads, -Q 20 sets the 
quality score (reads cut off if they fall below, but quality score is only used to break 
ties in the case of multiple overlaps), -C 70 sets the percentage of the read cut off if 
it falls below the quality score.

```> flash2 -M 200 -O -C 70 -Q 20 output_nodup_PE1.fastq output_nodup_PE2.fastq```

These were rather aggressive cleaning parameters but were chosen since many of the samples 
were sequenced quite deeply for their complexity. Flash2 output three files: 
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
combined: 23.61%. Flash2 also discarded a very small percentage (0.01%) due to low quality.
#### Trim Low Quality Ends
[Sickle](https://github.com/najoshi/sickle) was used to trim both the overlapped and 
non-overlapped reads output by flash2. Reads produced from most modern sequencing technologies 
have progressively lower quality approaching the 3'-end and sometimes the 5'-end as well. 
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
got run separately. First the paired end reads were processed using the sickle pe command:

```> sickle pe -n --length-threshold 75 --qual-threshold 20 --qual-type sanger -f 
out.notCombined_1.fastq -r out.notCombined_2.fastq -o cleaned_PE1.fastq -p cleaned_PE2.fastq 
-s cleaned_SE1.fastq```

The output is three files: **cleaned_PE1.fastq**, **cleaned_PE2.fastq**, and **cleaned_SE1.fastq**. 
There are two files for the forward and reverse reads and one for single reads. 
The single reads were created when the quality of one read in the pair was too low causing 
it to be discarded. The screen output is:

```PE forward file: out.notCombined_1.fastq 
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
removed from each end of the input reads.

Sickle se was run on the file of single overlapped reads output by flash2. It was run with
the same parameters as sickle pe:

```> sickle se -n --length-threshold 75 --qual-threshold 20 --qual-type sanger --fastq-file
 out.extendedFrags.fastq --output-file cleaned_SE2.fastq```
 
It output the file: **cleaned_SE2.fastq**, as well as the screen output:

```SE input file: out.extendedFrags.fastq

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

### Metagenomic Assembly
The cleaned reads were assembled into longer contiguous sequences (contigs) using the 
assembler [SPAdes version 3.9.0](http://cab.spbu.ru/software/spades/). It was run in meta
mode with the default k-mer sizes of 21, 33, and 55.

```> spades.py --meta -1 cleaned_PE1.fastq -2 cleaned_PE2.fastq -s cleaned_SE.fastq -o
spades_output/```

This allowed it to use up to 40 threads and 400Gb of memory and it took ~5.5 hours to run.
The output is a directory called **spades_output**. This includes a variety of files, the 
two most important of which are: **contigs.fasta** and **spades.log**. The file 
**spades.log** gives information about the run including how long the run took and whether
the run failed or finished. The commands:

```head -n 2 contigs.fasta 
>NODE_1_length_648482_cov_282.424
CCGAAGCCCGCCGATGCGGGCTTCGGCCGTTCCGGGCCCGCCTTTTCGGCGGGCCTTTGC
```
and:
```tail -n 2 contigs.fasta 
>NODE_3179_length_56_cov_607
CGGGGTGTGTAGGGCGAATAACGCCATGGGCGTTATCCGCCGATGTCGCGGATGAG
```

can be used to get a quick idea of how well the assembly worked. The reads from Community 4a
assembled into 3179 total contigs, the longest of which was 648482bp. The shortest contigs 
were only 56bp and not useful for further analysis. Only those contigs longer than 500bp
were kept for further analysis.




