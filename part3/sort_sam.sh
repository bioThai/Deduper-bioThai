#!/usr/bin/bash



# call samtools to sort sam file by 

input_sam_path="../part1/"
input_sam_file="test_input.sam"
output_sam="sorted_${input_sam_file}"

/usr/bin/time -v samtools sort -M -o $output_sam

exit