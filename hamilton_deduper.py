#!/usr/bin/env python

import argparse
import re

def get_args():
    parser = argparse.ArgumentParser(description="Deduplicating PCR duplicate sequences")
    parser.add_argument('-u', '--umi', type = str, help = 'designates file/path to file containing the list of UMIs', required = True)
    parser.add_argument('-f', '--file', type = str, help = 'designates absolute file path to sorted sam file input', required = True)
    parser.add_argument('-o', '--outfile', type = str, help = 'designates absolute file path to sorted sam file output', required = True)
    #parser.add_argument('-h', '--help', type = str, help = '_____ input help message here')
    return parser.parse_args()
args = get_args()

umi_file = args.umi
sam_file = args.file
output_file = args.outfile

def parse_UMI(UMI_file:str) -> list:
    '''Takes in a file of known UMIs where each UMI is its own line of the file
    and outputs a set of the known UMIs used'''
    known_UMIs = []
    lines = 0
    with open(UMI_file, "r") as fumi:
        for line in fumi:
            line = line.strip("\n")
            known_UMIs.append(line)
            lines += 1
        return known_UMIs
        


def parse_SAM_file(line:str) -> list:
    '''Takes in the sequence alignment information from a SAM file and returns a list containing 
the sequence UMI, Chromosome, Position, CIGAR string, and Strand'''
    #parse individual samfile line
    line = line.split("\t")
    #extract chromosome, position, cigar and bitwise flag from the samfile line
    Chrom = line[2]
    Pos = line[3]
    CIGAR = line[5]
    ID_line = line[0].split(":")
    #extract UMI
    UMI = ID_line[-1]
    flag = int(line[1])
    #check for strandedness within bitwise flag
    if ((flag & 16) == 16):
        Strand = "-"
    else: 
        Strand = "+"
    dup_info = [UMI, Chrom, Pos, CIGAR, Strand]        
    return dup_info


def parse_CIGAR(dup_info: list) -> list:
    '''Given an input list of UMI, chromosome, position, CIGAR string, and strand (+ or -) 
    for each SAM alignment, will determine the 5' start position of the read and return list:
        UMI, chromosome, 5' position, and strand'''
    
    CIGAR = dup_info[3]
    final_pos = int(dup_info[2])
    strand = dup_info[-1]

    # Parse CIGAR string to get operations and their lengths
    cig = re.findall(r"(\d+)([MIDNSHP=X])", CIGAR)

    if strand == "+":
        # Adjust final_pos by adding the length of any initial soft clipping (S)
        if cig[0][1] == "S":
            final_pos -= int(cig[0][0])
    
    elif strand == "-":
        # On the negative strand, calculate the position by adding all lengths
        # of matches (M), deletions (D), and skipped regions (N).
        for length, op in cig:
            if op in "MDN":
                final_pos += int(length)
        # Account for trailing soft clipping (S) if present
        if cig[-1][1] == "S":
            final_pos += int(cig[-1][0])

    # Return updated information with the computed 5' position
    return [dup_info[0], dup_info[1], final_pos, strand]


#create list of known UMIs from file
known_umis = parse_UMI(umi_file)

#sort for pcr duplicates
with open(sam_file, "r") as infile:
    with open(output_file, "w") as outfile:
        current_dupes = {}
        chr_count = 0

        for line in infile:
            line = line.strip()
            
            # Add header lines directly to the output file
            if line.startswith("@"):
                outfile.write(f"{line}\n")
                continue
            
            # Parse read information
            ID_info = parse_SAM_file(line)
            UMI = ID_info[0]
            chrom = ID_info[1]

            # If this is a new chromosome, reset current_dupes and set chr_count
            if chrom != chr_count:
                chr_count = chrom
                current_dupes = {}
                # Reset duplicates for the new chromosome

            # Only check position of reads with known UMIs
            if UMI in known_umis:
                adj_pos_ID_info = parse_CIGAR(ID_info)
                position = adj_pos_ID_info[2]
                strand = adj_pos_ID_info[3]
                
                # Create a unique identifier for each read
                read_ID = (UMI, position, strand)
                
                # Write the read to output if it is not a duplicate
                if read_ID not in current_dupes:
                    current_dupes[read_ID] = line
                    outfile.write(f"{line}\n")




                

                






