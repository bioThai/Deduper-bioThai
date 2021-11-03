# READ_ME

## Reference Based PCR Duplicate Removal Tool

The deduplicator (deduper) Python script in this repository takes in SAM files of single-end RNA-seq reads and a text file of unique molecular identifiers (UMIs) to be used as a reference. The deduper script outputs summary statistics about the output deduped SAM file, such as:
- Number of header lines
- Number of reads with UMI errors (eg, UMIs not matching the ones in the given UMI file)
- Number of duplicate reads removed
- Number of unique reads (total)
- Number of unique reads (by chromosome/contig)

Additional features currently under construction:
- Paired-end functionality
- Ability to specify which PCR duplicate to keep based on quality score

### How to Use

To deduplicate a SAM file using SLURM on an HPC cluster:

1. Specify the following in the `deduper_wrapper.sh` script:
    - `input_sam_dir`: path to the directory that contains the SAM file to be deduplicated (include one forward slash `/` at the end of the path)
    - `input_sam_file`: the name of the SAM file to be deduplicated
    - `umi_file`: the path to and filename of the text file containing the list of reference UMIs (eg, `"../STL96.txt"`)
    - SBATCH directives for your HPC cluster (if different from what is provided in the script).

2. Run the `deduper_wrapper.sh` script as a SLURM batch job using the following command:

    `sbatch deduper_wrapper.sh`

__Note:__ Multiple SAM files can be deduplicated in a single run if each file is sorted via Samtools. This can be done by modifying the `deduper_wrapper.sh` script to include extra lines/commands for each file to be sorted. Then, all of these sorted SAM files can be inputted as arguments (with a space between each file name) for the Python script.
