#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --output=sort_sam%j.out
#SBATCH --error=sort_sam%j.err

conda activate bgmp_py39

input_sam_path="/projects/bgmp/shared/deduper/"
#"/projects/bgmp/tnguye14/bioinfo/Bi624/Deduper-bioThai/part1/"
input_sam_file="C1_SE_uniqAlign.sam"
#"test_input.sam"
output_sam="sorted_${input_sam_file}"

# call samtools to sort sam file by chromosome and left-based read pos, using 8 threads
/usr/bin/time -v samtools sort -M -@ 8 -o $output_sam ${input_sam_path}${input_sam_file}

# command line version:
# samtools sort -M -o sorted_test_input.sam /projects/bgmp/tnguye14/bioinfo/Bi624/Deduper-bioThai/part1/test_input.sam
exit