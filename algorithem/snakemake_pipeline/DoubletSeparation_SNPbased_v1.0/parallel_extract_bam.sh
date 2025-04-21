#!/bin/bash

# Parameter checking
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <readID_file> <input_bam> <output_dir> <threads>"
    exit 1
fi

readID_file=$1
input_bam=$2
output_dir=$3
threads=$4

# Create output directory
mkdir -p "$output_dir"

# Create a temporary folder under output_dir
temp_dir="$output_dir/temp_dir"
mkdir -p "$temp_dir"

# Define the prefix path of the split file
split_prefix="$temp_dir/readID_part"

# Split ReadID files into temp_dir
split -l 200000 "$readID_file" "$split_prefix"

# Check whether the split file is generated successfully
if ! ls "$temp_dir/readID_part"* 1>/dev/null 2>&1; then
    echo "Error: No split files generated in $temp_dir"
    exit 2
fi

# Define the processing function
process_part() {
    part_file=$1
    part_output="$temp_dir/$(basename "$part_file").bam"  
    echo "Processing $part_file -> $part_output"
    samtools view -N "$part_file" -b "$input_bam" -o "$part_output"
    if [ $? -ne 0 ]; then
        echo "Error processing $part_file"
    fi
}

export -f process_part
export input_bam
export temp_dir

# Process all split files in parallel and store the output BAM files in temp_dir
find "$temp_dir" -type f -name "readID_part*" | parallel -j "$threads" process_part

# Check if the BAM file is generated
if ! ls $temp_dir/readID_part*.bam 1>/dev/null 2>&1; then
    echo "Error: No BAM files generated in $temp_dir"
    exit 3
fi

# Merge all generated small BAM files
output_bam="$output_dir/$(basename "$readID_file" .txt)_separated_bam.bam"
separated_bam="$temp_dir/readID_part*.bam"
samtools merge -f $output_bam $separated_bam

# Clean up temporary files
rm -rf "$temp_dir"
echo "All extracted reads saved to $output_bam"
