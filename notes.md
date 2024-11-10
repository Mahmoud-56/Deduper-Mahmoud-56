
The goal of this assignment is to write a reference based PCR duplicate removal tool. That is, given a sorted sam file of uniquely mapped reads, remove all PCR duplicates (retain only a single copy of each read)


# Problem 
The problem we are trying to address here is having PCR duplicates in our data. Those duplicates arise because not every molecule is amplified equally during library prep, leading to false representation and interpretation of the data, especially for downstream analysis after alignment (It does not significantly impact alignment, but it does impact differential expression analysis for example). 

When looking at a SAM file, PCR duplicates should map to the same chromosome (RNAME: SAM col 3), 5' start of read (POS: SAM col 4), and be on the same strand (BITFLAG 16: SAM col 2). In addition, each read header will contain a UMI which indicates the origin of the read. If two reads mapped to the same region and have the same UMI, then they will be considered as duplicates in our analysis. An important consideration when comparing the mapping position between two reads is whether any of the reads has been soft-clipped. Soft-clipping happens at the start or end of the sequence and will change the position of the read. We can tell by looking at the CIGAR string letter 'S'. If a read has been soft clipped, we need to account for that by subtracting the number of clipped nucleotides from the value of the position. 

# Pseudocode 

### Input:
- Sorted SAM file
- List of UMIs

### Output 
- Deduplicated SAM file
- Duplicated reads (maybe)
# Pseudocode

```
Initialize an empty set called `read_pos` to track unique read positions (UMI, chromosome, adjusted start position).

Open the sorted SAM file (sorted by chromosome number using samtools or awk).

For each line in the SAM file:
    Extract the UMI from the read header.
    
    If UMI is not in the known UMI list:
        Discard the read or write it to a separate file.
        Continue to the next read.
    
    Extract the chromosome number and starting position (adjust for soft clipping based using CIGAR string).
    
    If we encounter a new chromosome number:
        Clear the `read_pos` set to start fresh for the new chromosome.
    
    Check if the read is on the reverse strand (using BITFLAG 16 from column 2):
        Adjust the start position based on the strand direction 
    
    Create a tuple (UMI, chromosome, corrected start position)

    If the tuple is already in read_pos:
        Discard the read (consider it a duplicate).
    Else:
        Add the tuple to read_pos.
        Write the read to the output file.

Close the SAM file and output file.

```

# High-level Function 

### 1. Find chromosome number
```
def find_chr_num(line: str) -> : int

'''takes in a the line containing the read info and outputs the number of chromosome where that read belongs'''
	return chr_num
	
Input:
NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	0	3	76814284	36	71M	*	0	0	TCCACCAC.... 

Output:
3
```

### 2. Find UMI 
```
def find_UMI(line: str) -> : str

'''takes in a the line containing the read info and outputs its UMI'''
	return UMI
	
Input:
NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	0	3	76814284	36	71M	*	0	0	TCCACCAC.... 

Output:
CTGTTCAC
```
### 3. Find Strand Orientation
```
def find_orientation(line: str) -> : str

'''takes in a the line containing the read info and outputs its orientation'''
	return orientation
	
Input:
NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	0	3	76814284	36	71M	*	0	0	TCCACCAC.... 

Output:
Plus

Input:
NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	16	3	76814284	36	71M	*	0	0	TCCACCAC.... 

Output:
Minus 
```

### Find Read Position 

```
def find_pos(line: str) -> : int

'''takes in a the line containing the read info and outputs its position'''
	return pos
	
Input:
NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	0	3	76814284	36	71M	*	0	0	TCCACCAC.... 

Output:
76814284

Input:
NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	0	3	76814286	36	2S69M	*	0	0	TCCACCAC.... 

Output:
76814284
```



### Part 3 

#### Output 
```
./Dedy.sh 
Number of duplicates removed: 4467362
Number of wrong UMIs found: 0
	Command being timed: "./AlMahmoud_Deduper.py -u STL96.txt -f C1_SE_uniqAlign_sorted.sam -o Dedup_output_SE"
	User time (seconds): 73.35
	System time (seconds): 5.06
	Percent of CPU this job got: 98%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:19.87
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 639372
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 162341
	Voluntary context switches: 1118
	Involuntary context switches: 143
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

---------------------------------------------------------------------------------------
using samtools to sort by chromosome number (default) before running python script 

```
samtools sort -o <output> <input> 
```

To get the counts of reads per chromosome, I used the following awk command:
```
awk '!/^@/ {print $3}' Dedup_output_SE | sort -V | uniq -c | awk '{print $2 "\t" $1}'
```


























