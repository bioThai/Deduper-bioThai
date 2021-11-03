#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=4
#SBATCH --time=10:00:00
#SBATCH --output=deduper_wrapper_%j.out
#SBATCH --error=deduper_wrapper_%j.err

conda activate bgmp_py39

input_sam_dir="/projects/bgmp/shared/deduper/"
input_sam_file="C1_SE_uniqAlign.sam"
sorted_sam_file="sorted_${input_sam_file}"
umi_file="../STL96.txt"

# call samtools to sort sam file by chromosome and left-based read pos, using 4 threads
# sorted SAM file will be created in current working directory
/usr/bin/time -v samtools sort -@ 4 -o $sorted_sam_file ${input_sam_dir}${input_sam_file}

# deduplicate sorted SAM file using a text file containing UMIs for reference
# output deduped SAM file will be in current working directory
/usr/bin/time -v python nguyen_deduper.py -f $sorted_sam_file -u $umi_file

#remove sorted SAM file from current working directory to save storage space
rm $sorted_sam_file

exit