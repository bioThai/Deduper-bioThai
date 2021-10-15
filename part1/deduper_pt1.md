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
- conda activate bgmp_py39 in the wrapper for Python script
- Import modules including pysam (the python equivalent of command line samtools).

- Get argparse inputs (filenames, etc.) and save them into variables.
- Create string variable to hold name of sorted version of input SAM file (eg, sorted_input.sam)
- Create string variable to hold name of deduplicated output SAM file (eg, input_deduped.sam)
- Create temp_read dictionary to temporarily hold one sam file line
    - # key: tuple containing read_umi, chromosome_name, starting_position
    - # value: list containing bitflag and rest of the values in the line
- Create a list to hold correct UMIs

- With open correct_UMIs.txt in read mode:
    - For line in file:
        - strip "\n" from line
        - Append line to correct_umis_list

- With open input.sam file in read mode:
    - Use pysam's sort function to sort input.sam file by chromosome name and read position ("-M" option) and write the sorted output into another file called sorted_input.sam (-o option). This way, the chromosomes and read positions are in order, so any duplicate reads are closer together in the file, causing less intensive memory usage when parsing file.

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
                    - starting_position = starting_position - bps_clipped
                
                - Create tuple containing read_umi, chromosome_name, starting_position, and assign it to "key" variable
                - If temp_read dictionary is empty:
                    - temp_read[key] = list(bitflag) + line_tokens[4:]
                - Else:
                    - If key is in temp_read dictionary:

                        Work on this more to develpoe dict that truly hold unique, non-duplicate read
                        - If bitflag of current read and bitflag of temp_read[key] (stored in temp_read[key][0])
                        
                        
                        next read have the same bit-16 flag states (either BOTH have bit 16 flipped, or BOTH don't have bit 16 flipped)

                
                    - at some point, check if bitwise flag of one read and bitwise flag of next read have the same bit-16 flag states (either BOTH have bit 16 flipped, or BOTH don't have bit 16 flipped):
                        - if they have the same bit-16 flag state, then they both have the same strandedness and possibly could be duplicates 
                
                - 
            

            - Else:
                - continue, read next line in sorted_input file
            

           
            
            
                    

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
    parser.add_argument("-f", nargs="+", help="specifies input SAM filename", type=str, required=True)
    parser.add_argument("-u", nargs="+", help="specifies input text file containing all possible UMIs", type=str, required=True)
    return parser.parse_args()

def is_valid_umi(str: umi_read) -> boolean:
    '''Checks if a given UMI is valid (error-free) by seeing if it matches any UMIs in a given dictionary '''


```

For this portion of the assignment, you should design your algorithm for single-end data, with 96 UMIs. UMI information will be in the QNAME, like so: ```NS500451:154:HWKTMBGXX:1:11101:15364:1139:GAACAGGT```. Discard any UMIs with errors (or error correct, if you're feeling ambitious).