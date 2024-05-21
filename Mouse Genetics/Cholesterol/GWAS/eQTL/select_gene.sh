#!/bin/bash

# Check if correct number of arguments are passed
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <file.tsv> <header> <output.txt>"
    exit 1
fi

file=$1
header=$2
output=$3

# Extract the column and skip the header
csvcut -t -c "$header" "$file" | tail -n +2 > "$output"