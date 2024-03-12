#!/bin/bash

bedFile=$1

bed12File=${bedFile%.*}.bed12

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6, $2, $3, "0", "1", $3 - $2, "0"}' "$bedFile" > "$bed12File"

echo "Converted bed6 to bed12: $bedFile >>> $bed12File"
