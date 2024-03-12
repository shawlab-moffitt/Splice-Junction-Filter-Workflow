#!/bin/bash

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 input.bed|bed_list.txt num_nucleotides"
    exit 1
fi

input="$1"
num_nucleotides="$2"

process_bed() {
    local input_file="$1"
    while IFS=$'\t' read -r chrom start end rest; do
        new_start=$((start - num_nucleotides))
        new_start=$((new_start < 0 ? 0 : new_start))
        
        new_end=$((end + num_nucleotides))
        
        echo -e "$chrom\t$new_start\t$start\t$rest" >> "$output_before"
        echo -e "$chrom\t$end\t$new_end\t$rest" >> "$output_after"

    done < "$input_file"
    
    echo "Processed: $input_file"
}

if [ -f "$input" ]; then
    num_columns=$(awk -F'\t' '{print NF; exit}' "$input")
    if [ "$num_columns" -eq 1 ]; then
        while IFS= read -r bed_file; do
            FileNoExt=${bed_file%.*}
            output_before="${FileNoExt}_${num_nucleotides}LeftRegion.bed6"
            output_after="${FileNoExt}_${num_nucleotides}RightRegion.bed6"
            process_bed "$bed_file"
        done < "$input"
    elif [ "$num_columns" -eq 6 ]; then
        FileNoExt=${input%.*}
        output_before="${FileNoExt}_${num_nucleotides}LeftRegion.bed6"
        output_after="${FileNoExt}_${num_nucleotides}RightRegion.bed6"
        process_bed "$input"
    else
        echo "Input format not recognized."
        exit 1
    fi
else
    echo "Input file not found."
    exit 1
fi


