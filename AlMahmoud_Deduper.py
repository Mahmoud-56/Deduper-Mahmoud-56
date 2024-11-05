#!/usr/bin/env python

import argparse
import re

def get_args():
    parser = argparse.ArgumentParser(description="Python program for reference based PCR duplicate removal.")
    parser.add_argument("-f", "--input", help="file path to SAM file.", type=str, required=True)
    parser.add_argument("-o", "--outfile", help="file path to deduplicated SAM file.", type=str, required=True)
    parser.add_argument("-u", "--umi", help="file containing a list of UMI  barcodes.", type=str, required=True)
    return parser.parse_args()



def get_UMI_list(UMIs: str) -> set:
    '''Takes in a file containing a list of one barcode per line, it returns a set of those barcodes'''
    umis_set = set()
    with open(UMIs, 'r') as barcode_file:
        while True:
            umis = barcode_file.readline().strip()  # remove white space
            if not umis:  # if the line is empty, break the loop
                break
            if set(umis).issubset("AaGgTtCc"):  # check if it is a valid barcode
                umis_set.add(umis)  # add the cleaned-up UMI to the set
            else:  # if the barcode has an invalid base, it is not added to the set
                print(f"Invalid UMI: {umis}")
    return umis_set


def parse_sam_line(line: str):
    '''Extracts UMI, chromosome, strand, CIGAR, and position from a SAM line'''
    if line.startswith('@') or not line.strip():  # Check for header lines or empty lines:
        return None
    
    fields = line.split('\t')

    if len(fields) < 11:  #Missing column is SAM file
        
        return None
    
    read_name = fields[0]
    chromosome = fields[2]
    pos = int(fields[3])
    cigar = fields[5]
    flag = int(fields[1])

    # Determine strand from bitwise FLAG column 2
    if (flag & 16) == 16:
        strand = '-' 
    else:
        strand = '+'

    # Extract UMI from read name (UMI is in the last section of read name after a colon such as readname:UMI)
    umi = read_name.split(':')[-1]

    return umi, chromosome, strand, cigar, pos

def adjust_5prime_start_position(cigar: str, pos: int, strand: str) -> int:
    '''Adjusts the 5' starting position based on the CIGAR string and strand direction'''
    if strand == '+':
        soft_clip_match = re.match(r"(\d+)S.+", cigar) #find if there is an S in the CIGAR 
        if soft_clip_match:
            pos -= int(soft_clip_match.group(1))  # adjust position for soft clipping (group 1 captures the first match of soft clipping -> only left side)
    else: # - strand
        # Sum M, D, N, and consider any soft clipping on the right side
        matches = re.findall(r"(\d+)[MDN]", cigar)  # find all matches of M, D, and N ignoring soft clipping for now 
        #add each value to the pos number
        pos += sum(map(int, matches)) # map applies the int function to every match avoiding using a loop

        right_soft_clip = re.search(r"(\d+)S$", cigar) #look for soft clipping on the right side only 
        if right_soft_clip:
            pos += int(right_soft_clip.group(1))  # adjust right side soft clipping
    return pos

 


def remove_duplicates(input_sam: str, output_sam: str, known_umis: set):
    '''Removes PCR duplicates from a sorted SAM file based on UMI, strand, and adjusted start position'''
    current_chromosome = None
    unique_reads = set()  # Track unique (UMI, strand, adjusted position) combinations
    duplicate_count = 0
    wrong_umis_count = 0

    with open(input_sam, 'r') as input, open(output_sam, 'w') as output:
        for line in input:
            if line.startswith('@'):
                output.write(line)
                continue

            IDs = parse_sam_line(line)
            if not IDs:  # Check if IDs is None, meaning the line was invalid or a header
                continue

            umi, chromosome, strand, cigar, pos = IDs

            # Check if UMI is valid
            if umi not in known_umis:
                wrong_umis_count += 1    
                continue

            # Clear unique reads set if new chromosome is encountered 
            if chromosome != current_chromosome:
                current_chromosome = chromosome
                unique_reads.clear()

            # Adjust position based on CIGAR string
            adjusted_pos = adjust_5prime_start_position(cigar, pos, strand)

            # Identify read by UMI, strand, and adjusted position
            read_identifier = (umi, strand, adjusted_pos)

            # Write if unique; otherwise skip duplicates
            if read_identifier not in unique_reads:
                unique_reads.add(read_identifier)
                output.write(line)
            else:
                duplicate_count += 1

    return duplicate_count, wrong_umis_count



def main():
    args = get_args()
    known_umis = get_UMI_list(args.umi)
    duplicate_count, wrong_umis_count = remove_duplicates(args.input, args.outfile, known_umis)
    print(f"Number of duplicates removed: {duplicate_count}")
    print(f"Number of wrong UMIs found: {wrong_umis_count}")

if __name__ == "__main__":
    main()

