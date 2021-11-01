#!/bin/python
import argparse
import re

def get_args():
    '''Defines/sets possible command line arguments for script'''
    parser = argparse.ArgumentParser("A program to parse SAM file and filter out reads that are PCR duplicates or that have UMIs which don't match a given list of possible UMIs.")
    parser.add_argument("-f", "--file", nargs="+", help="Specifies input SAM filename(s). SAM files must be unzipped and must already be sorted by chromosome name (RNAME) and 1-based leftmost mapping position (POS).", type=str, required=True)
    parser.add_argument("-p", "--paired", help="Specifies if input SAM file holds paired-end alignment data. If paired-end, then specify the Boolean value 'true' (without quotations) as argument when calling this script. If not paired-end, do not include -p option.", type=bool)
    parser.add_argument("-u", "--umi", help="Specify the text file(s) containing the list of UMIs to compare input SAM file UMIs to. Each UMI must be on a different line and not be delimited by any characters. If this option is unset/unspecified, then randomers will be used for UMIs instead of an UMI text file.", type=str)
    
    #print error messages?
    #parser.add_argument("-h", "--help", help="")
    return parser.parse_args()

def get_umi_list(umi_file: str) -> list:
    '''Takes a text file of UMIs and reads them into a list. Returns the list of UMIs.'''
    #local variables
    umi_list: list = list()

    #read text file of UMIs and put each UMI into a list
    with open(umi_file, "r") as fh:
        for line in fh:
            line = line.strip("\n")
            umi_list.append(line)
    return umi_list

def is_positive_strand(bitflag: int) -> bool:
    '''Takes a bitwise flag from a SAM file line and returns whether the read associated with that bitwise flag belongs to a positive or negative strand.'''
    #local variables
    is_pos_strand: bool = True

    if ((bitflag & 16) == 16):
        is_pos_strand = False
    else:
        is_pos_strand = True   
    return is_pos_strand     

def get_adjusted_start_pos(sam_start_pos: int, cigar_string: str, is_pos_strand: bool) -> int:
    '''Takes the left-most start position (according to SAM file), CIGAR string, and strandedness of a SAM file read and calculates the "true" start position of the read.'''
    # local variables
    adjusted_start_pos: int = 0
    left_clipped: int = 0    #holds number of base pairs soft-clipped from left-most position of read
    cigar_search_dict: dict = {"M": 0, "D": 0, "S": 0, "N": 0}
        #holds number of bases with Ms (matches or mismatches), Ds (deletions), right-most Ss (right-most soft-clipped bases), and Ns (intron sequences, which are present in reference genome but not in the RNA-seq reads) from cigar string
        #key: "M", "D", "S", "N" (Ms, Ds, right-most Ss, and Ns affect the "true" start position of a negative strand)
        #value: total number of bases
    temp_num: int = 0

    if is_pos_strand == True:
        # extract number of basepairs soft-clipped from beginning of read position
        left_clipped_search_results_list = re.findall("^([0-9]*)(S)", cigar_string)
        if len(left_clipped_search_results_list) != 0:
            left_clipped = int(left_clipped_search_results_list[0][0])
            adjusted_start_pos = sam_start_pos - left_clipped
        else:
            adjusted_start_pos = sam_start_pos
    #else if sequencing read currently being read is from a negative strand
    else:
        # use regex with f-string to find number of Ms, Ds, Ss, and Ns from cigar string
        for letter in cigar_search_dict:
            search_results_list = re.findall(f"([0-9]*)({letter})", cigar_string)
            if len(search_results_list) != 0:
                temp_num = 0
                for regex_match in search_results_list:
                    temp_num += int(regex_match[0])
                cigar_search_dict[letter] = temp_num
                if letter == "S":
                    # if read is from negative strand, then leftmost soft-clipping doesn't affect true read start position, 
                    # so you only want the right-most soft-clipping (total soft clippng - left-sided soft clipping)
                    left_clipped_search_results_list = re.findall("^([0-9]*)(S)", cigar_string)
                    if len(left_clipped_search_results_list) != 0:
                        left_clipped = int(left_clipped_search_results_list[0][0])
                        cigar_search_dict[letter] = temp_num - left_clipped
            else:
                cigar_search_dict[letter] = 0
        adjusted_start_pos = sam_start_pos + sum(cigar_search_dict.values()) - 1
        
    return adjusted_start_pos
    

def dedup_sam(input_sam_file: str, umi_list: list, is_paired_end: bool):
    '''Takes a reference list of valid UMIs and a sorted input SAM file, and deduplicates the SAM file.'''
    #local variables
    output_sam_file: str = input_sam_file.split(".sam")[0] + "_deduped.sam"  
        #create an output sam file name based on the input sam file name
    temp_reads_dict: dict = {}  
        # holds all non-duplicated reads from a specific chromosome
        # key: tuple containing read_umi, bit16_flag (true or false), chromosome_name, starting_position
        # value: list containing all items from line read from SAM file
    current_chrom_in_dict: str = ""     #holds value of which chromosome's non-duplicate reads are currently being stored in temp_reads dictionary
    line_tokens: list = []              #holds each column from an aligned sequencing line in SAM file as a list element
    qname: str = ""
    bitflag: int = 0
    chrom_name: str = ""
    sam_starting_pos: int = 0
    cigar_string: str = ""
    read_umi: str = ""
    is_pos_strand: bool = True #if current read is from a positive strand (if false, read is from negative strand)
    adjusted_start_pos: int = 0
    key: tuple = ()
    output_line_tokens: list = []
    output_read: str = ""
    num_header_lines: int = 0           #holds count of number of header lines in input SAM file
    num_umi_error_reads: int = 0
    num_duplicate_reads: int = 0
    num_unique_reads_dict: dict = {}
        #holds the number of unique, non-duplicate reads for each chromosome/contig
        #key: chromosome/contig name
        #value: number of unique reads

    with open(input_sam_file, "r") as input_sam_fh, open(output_sam_file, "w") as output_sam_fh:
        for line in input_sam_fh:
            line = line.strip("\n")
            #if line is a header line in sam file (starts with "@")
            if line.startswith("@"):
                num_header_lines += 1
                output_sam_fh.write(line + "\n")
            else:
                line_tokens = line.split()
                qname = line_tokens[0]
                bitflag = int(line_tokens[1])
                chrom_name = line_tokens[2]
                sam_starting_pos = int(line_tokens[3])
                cigar_string = line_tokens[5]

                read_umi = qname.split(":")[-1]

                if read_umi not in umi_list:
                    num_umi_error_reads += 1
                    continue
                else:
                    #check if current read is for a negative or positive strand
                    is_pos_strand = is_positive_strand(bitflag)
                    #get adjusted start position for current read based on its strandedness and corresponding cigar string
                    adjusted_start_pos = get_adjusted_start_pos(sam_starting_pos, cigar_string, is_pos_strand)
                    #create tuple containing 4 properties of a read (if 2 reads have the exact same values for these 4 properties, then they are PCR duplicates)
                    key = (read_umi, is_pos_strand, chrom_name, adjusted_start_pos)

                    if len(temp_reads_dict) == 0:
                        current_chrom_in_dict = chrom_name
                        num_unique_reads_dict[current_chrom_in_dict] = 0

                    #if current read is a PCR duplicate
                    if key in temp_reads_dict:
                        num_duplicate_reads += 1
                        continue
                    else:
                        if chrom_name == current_chrom_in_dict:
                            temp_reads_dict[key] = line_tokens
                            num_unique_reads_dict[current_chrom_in_dict] += 1
                        else:
                            for non_duplicated_read_key in temp_reads_dict:
                                #Get the corresponding value (a list) out of dict and save it in output_line_tokens
                                output_line_tokens = temp_reads_dict[non_duplicated_read_key]
                                #Join all values in output_line_tokens list into a single string, using a "\t" as delimiter
                                output_read = "\t".join(output_line_tokens)
                                output_sam_fh.write(output_read + "\n")
                            #clear the temp reads dict before adding reads from new chromosome
                            temp_reads_dict.clear()
                            #add current read (first read from new chromosome) to dict
                            temp_reads_dict[key] = line_tokens
                            current_chrom_in_dict = chrom_name
                            #need to initialize new entry in unique reads dict with value of 0 for new chromosome encountered
                            num_unique_reads_dict[current_chrom_in_dict] = 0
                            num_unique_reads_dict[current_chrom_in_dict] += 1

        #Once all the non-duplicated reads from the last chromosome have been saved into temp_reads dict, write these last reads into output file:
        for non_duplicated_read_key in temp_reads_dict:
            #Get the corresponding value (a list) out of dict and save it in output_line_tokens
            output_line_tokens = temp_reads_dict[non_duplicated_read_key]
            #Join all values in output_line_tokens list into a single string, using a "\t" as delimiter
            output_read = "\t".join(output_line_tokens)
            output_sam_fh.write(output_read + "\n")

    #print summary stats from de-duplication
    print("Input file:", input_sam_file)
    print("Output file:", output_sam_file)
    print("Number of header lines:", num_header_lines)
    print("Number of reads with UMI errors:", num_umi_error_reads)
    print("Number of duplicate reads:", num_duplicate_reads)
    print("Number of unique reads (total):", sum(num_unique_reads_dict.values()))
    print("\nNumber of unique reads (by chromosome/contig):")
    for chromosome in sorted(num_unique_reads_dict.keys()):
        print(chromosome, num_unique_reads_dict[chromosome], sep="\t")


def main():
    '''Main function, drives the order of execution for script'''
    #local variables
    args = get_args()
    input_sam_file_list: list = args.file
    input_umi_file: str = args.umi
    is_paired_end: bool = args.paired
    umi_list: list = []

    #call function to read input_umi_file and save UMIs into list
    umi_list = get_umi_list(input_umi_file)
    #print(len(umi_list), umi_list)
    
    for input_sam_file in input_sam_file_list:
        dedup_sam(input_sam_file, umi_list, is_paired_end)
        

if __name__ == "__main__":
    main()
