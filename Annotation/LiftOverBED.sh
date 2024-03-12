#!/bin/bash



inFile=$1
inGenome=$2
LiftOverRef=$3
outPrefix=$4
outGenome=$5




tempBedName="${outPrefix}_TEMPlift.bed"
mappedBedName="${outPrefix}_MappedLift.bed6"
unmappedBedName="${outPrefix}_unMappedLift.bed6"
tempJxnListFile="${outPrefix}_TEMPLiftJxnList.txt"
JxnListFile="${outPrefix}_LiftJxnList.txt"
echo -e "Junctions_${inGenome}Build\tJunctions_${outGenome}LiftBuild" >> "${tempJxnListFile}"

process_LiftOver() {
	local input_file="$1"
	local newBed_file="$2"
	local JxnList_file="$3"
	while IFS=$'\t' read -r chrom start end name score strand
	do
		newName="${chrom}:${start}-${end}:${strand}"
		echo -e "${chrom}\t${start}\t${end}\t${newName}\t${score}\t${strand}" >> "${newBed_file}"
		echo -e "${name}\t${newName}" >> "${JxnList_file}"
	done < "$input_file"
}

liftOver ${inFile} ${LiftOverRef} ${tempBedName} ${unmappedBedName}

process_LiftOver ${tempBedName} ${mappedBedName} ${tempJxnListFile}
rm ${tempBedName}

if [[ $(wc -l <${unmappedBedName}) -gt 0 ]]
then
	sed '/^#/d' ${unmappedBedName} > unmappedBed_temp.txt
	cut -f4 unmappedBed_temp.txt > unmappedBed_temp2.txt
	sed -i "s/$/\tNA/" unmappedBed_temp2.txt
	cat ${tempJxnListFile} unmappedBed_temp2.txt > ${JxnListFile}
	rm ${tempJxnListFile} unmappedBed_temp.txt unmappedBed_temp2.txt
else 
	mv ${tempJxnListFile} ${JxnListFile}
fi

