# Deduper Lab Notebook

## 11/9/2024

I have had literally no time to work on this assignment so here we go! Deduper. Here is my pseudocode with attempts to implement the feedback given to me by peers:

```bash
def process_sam_file(file_path: str):
    Output: known_umis = set()  # list of known UMIs from STL96.txt
    
    
    with open(file_path, 'r') as sam_file, open('output.sam', 'w') as output_file:
        for line in sam_file:
            if line.startswith('@'):  # Header lines
                write_read_to_output(line, output_file)
                continue

            qname, chrom, pos, cigar, seq, qual = parse_sam_line(line)
            umi = extract_umi(qname)

            if not is_valid_umi(umi, known_umis):
                continue  # Skip invalid UMI reads

            read_hash = generate_read_hash(chrom, pos, umi)
            read_cigar_string = parse_CIGAR(read)
            
            if read_hash not in seen_reads:
                write_read_to_output(line, output_file)
                seen_reads.add(read_hash)
		write output to file

def parse_SAM_file(header: str) -> list:
"Takes in the sequence alignment information from a SAM file and returns a list containing 
the sequence UMI, Chromosome, Position, CIGAR string, and Position"
Input: '''NS500451:154:HWKTMBGXX:1:11101:11995:1145:AGCTACCA	0	2	76710746	36	71M	*	0	0	
TGCATAACTCGTGCTGGTTTCCTCCTTTGTGGGGACGTGATAGGTCGAGTACCTGAAGTCTCTTCTTCTGT	
6<EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	
MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU'''
	return ["AGCTACCA", 2, 76710746, 71M, "+"]

def generate_read_hash(chrom: str, pos: int, umi: str, strand: str) -> str:
    return f'{chrom(col3)}_{pos(col4)}_{umi(col1)}_{strand(col2, FLAG)}'

def parse_CIGAR(cigar: str) -> str:
		return start position of read

```

for bitwise flag:

bit 16 — if strand has to be rev. complemented it is from the - strand. So if 16 is true, minus strand is true. if reads on different strands, not PCR duplicates. 

## 11/10/24 and 11/11/24

Continued working on script

Even though I understand the concept behind fwd/rev position adjustment it is still difficult to accomplish. Even though I have something that works it’s slow, I need to keep working on it!

## 11/12/24

[carlyham@n0349 Deduper-carlyham]$ ./hamilton_deduper.py -f /projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam -u STL96.txt -o ../deduper/C1_SE_uniqAlign_output.sam
Unique Reads = 18186360

This is definitely not correct, I need to keep working on my script. currently the total # of reads is being listed as unique reads.

Note: IGV-fun environment has samtools in it . used it to sort files. 

[carlyham@login2 Deduper-carlyham]$ srun -A bgmp -p bgmp --time=0-5 ./hamilton_deduper.py -f ../deduper/sorted_by_chr.sam -u STL96.txt -o ../deduper/after_sort_attempt.txt
srun: job 23252556 queued and waiting for resources
srun: job 23252556 has been allocated resources
Unique Reads = 13719048

After sorting the file I am getting an answer that makes more sense. 

Highlights from slurm script of C1_SE_uniqAlign_output.sam:

```bash
Unique Reads = 13719048
	Command being timed: "./hamilton_deduper.py -f ../deduper/sorted_by_chr.sam -u STL96.txt -o ../deduper/after_sort_attempt.txt"
	Percent of CPU this job got: 98%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:06.82
	Maximum resident set size (kbytes): 1416488
	Exit status: 0
```

grep for header lines:  grep "^@" after_sort_attempt.txt | wc -l
65