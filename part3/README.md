# READ_ME

## Reference Based PCR Duplicate Removal Tool

The deduplicator script in this repository takes in SAM files of single-end RNA-seq reads and a text file of unique molecular identifiers (UMIs) to be used as a reference. 

Additional features currently under construction:
- Paired-end functionality
- Ability to specify which PCR duplicate to keep based on quality score

To deduplicate a SAM file using SLURM on an HPC cluster:

1. Specify the following in the `deduper_wrapper.sh` script:
    - `input_sam_dir`: path to the directory that contains the SAM file to be deduplicated (include one forward slash `/` at the end of the path)
    - `input_sam_file`: the name of the SAM file to be deduplicated
    - `umi_file`: the path to and filename of the text file containing the list of reference UMIs (eg, `"../STL96.txt"`)
    - SBATCH directives for your HPC cluster (if different from what is provided in the script).

2. Run the `deduper_wrapper.sh` script as a SLURM batch job using the following command:

    `sbatch deduper_wrapper.sh`
