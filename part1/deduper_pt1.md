# Deduper

## Part 1

Write up a strategy for writing a Reference Based PCR Duplicate Removal tool. That is, given a sam file of uniquely mapped reads, remove all PCR duplicates (retain only a single copy of each read). Develop a strategy that avoids loading everything into memory. You should not write any code for this portion of the assignment. Be sure to:

- Define the problem:
    - During the library prep process of RNA-Seq experiments, the cDNA libraries must be amplified by PCR so that detectable/sequenceable quantities are produced. However, the process of PCR amplification introduces noise and uncertainty into the data. This makes it harder to tell if duplicate sequencing reads are truly representative of gene expression levels in the original samples (i.e. if observed duplicates are due to natural biological duplicates indicative of higher expression levels), or if they're not truly representative of expression levels (i.e. if observed duplicates are  actually artificially-introduced PCR duplicates, which lead to inaccurate representations of gene expression levels).

- Write examples:
    - Include a properly formated input sam file: ```Deduper-bioThai/part1/test_input.sam```
    - Include a properly formated expected output sam file: ```Deduper-bioThai/part1/test_output.sam```

- Develop your algorithm using pseudocode:

```
- Import modules, including pysam (the python equivalent of command line samtools).

- Get argparse inputs (filenames, etc.) and save them into variables.

- Create string variable to hold name of sorted version of input SAM file (eg, sorted_input.sam)
- Create string variable to hold name of deduplicated output SAM file (eg, input_deduped.sam)
- Create temp_read dictionary to temporarily hold all non-duplicated reads from a specific chromosome
    - # key: tuple containing read_umi, bit16_flag (true or false), chromosome_name, starting_position
    - # value: list containing all items from line read from SAM file
- Create a list to hold correct UMIs

- With open correct_UMIs.txt in read mode:
    - For line in file:
        - strip "\n" from line
        - Append line to correct_umis_list

- With open input.sam file in read mode:
    - Use pysam's sort function to sort input.sam file by chromosome name and read position ("-M" option), and write the sorted output into another file called sorted_input.sam (-o option). 
    - This way, the chromosomes and read positions are in order, so any duplicate reads are closer together in the file, causing less intensive memory usage when parsing file.

- With open sorted_input.sam in read mode, with open input_deduped.sam in write mode:
    - For line in sorted_input filehandler:
        - strip "\n" from line
        - If line starts with "@" (line is a header line):
            - write line to output file (input_deduped.sam)
        - Else:  
            - split line by whitespace, save resulting list of line tokens into a variable called line_tokens 
            - save first line token into qname string variable
            - save second line token into bitflag int variable
            - save third line token into chromosome_name string variable
            - save fourth line token into starting_position int variable
            - save sixth line token into cigar_string variable

            - Split qname string by ":", and save last token into variable called read_umi
            - If read_umi is in correct_umis_list:
                - If "S" is in cigar_string (soft clipping occurred):
                    - Use regex to extract how many basepairs were soft-clipped, and save this in a variable called bps_clipped
                    - adjusted_starting_position = starting_position - bps_clipped
                - Else:
                    - adjusted_starting_position = starting_position
                - Check if bit 16 in bitflag is in "true" (1) or "false" (0) state:
                    - If bit 16 is true:
                        assign "true" to bit16_flag variable.
                    - Else:
                        assign "false" to bit16_flag.

                - Create tuple containing read_umi, bit16_flag, chromosome_name, adjusted_starting_position, and assign it to "key" variable. 
                    - Note: If bitwise flag of one read and bitwise flag of another read have the same bit-16 flag states (either BOTH have bit16_flag == TRUE, or BOTH have bit16_flag == FALSE), then they both have the same strandedness and possibly could be duplicates.

                - If temp_read dictionary is empty:
                    - Assign chromosome_name to current_chrom_in_dict variable (holds value of which chromosome's non-duplicate reads are currently being stored in temp_reads dictionary).

                - If key is in temp_read dictionary (current read is a pcr duplicate):
                    - continue to next iteration in for loop, read next line in sorted_input file
                - Else:
                    - If chromosome_name == current_chrom_in_dict:
                        - add read to dict (temp_read[key] = line_tokens)
                    - Else: 
                        - For key in temp_read dict:
                            - Pop the corresponding value (a list) out of dict and save it in an output_line_tokens list variable (output_line_tokens = pop(key))
                                - Join all values in output_line_tokens list into a single string, using a "\t" as delimiter. Assign this string to a variable.
                                - Write this string into the output sam file (input_deduped.sam).          
            - Else:
                - continue to next iteration in for loop, read next line in sorted_input file

```

- Determine high level functions:
    - Description
    - Function headers
    - Test examples (for individual functions)
    - Return statement
    
```
def main():
    '''Drives order of execution for script.'''
    # Calls get_args() to retrieve input arguments
    # calls functions to run rest of deduper script
    # No return statement

def get_args():
    '''Defines/sets possible command line arguments for script'''
    parser = argparse.ArgumentParser("A program to parse SAM file and filter out reads that are PCR duplicates or that have UMIs which don't match a given list of possible UMIs")
    parser.add_argument("-f", nargs="+", help="specifies input SAM filename(s)", type=str, required=True)
    parser.add_argument("-u", nargs="+", help="specifies input text file(s) containing all possible UMIs", type=str, required=True)
    return parser.parse_args()

def get_adjusted_start_pos(sam_start_pos: int, cigar_string: str) -> int:
    '''Takes a read start position and checks if corresponding cigar string has "S" in it(indicating soft clipping). If cigar string has S, then use regular expressions to extract the number of base pairs that were soft-clipped, and subtract this from the original start position value. Returns the adjusted start position value.'''

    Returns the adjusted start position value as integer
    
    # Examples:
    # get_adjusted_start_pos(100, 3S14M) returns 97
    # get_adjusted_start_pos(100, 14M20N30M) returns 100

```

For this portion of the assignment, you should design your algorithm for single-end data, with 96 UMIs. UMI information will be in the QNAME, like so: ```NS500451:154:HWKTMBGXX:1:11101:15364:1139:GAACAGGT```. Discard any UMIs with errors (or error correct, if you're feeling ambitious).