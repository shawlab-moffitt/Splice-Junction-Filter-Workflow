#!bin/bash

#condaEvn=$1
bedFile=$1
genome=$2
outDir=$3

mkdir -p ${outDir}

bedCols=$(awk -F'\t' '{print NF; exit}' ${bedFile})
if [ ${bedCols} -gt 1 ]
then
	fileName=${bedFile%.*}
	FileNameNoExt=$(basename ${fileName})
	newFileName="${FileNameNoExt}_HomerAnno.txt"
	annotatePeaks.pl ${bedFile} ${genome} > ${outDir}/${newFileName}
elif [ ${bedCols} -eq 1 ]
then
	while read i
	do
		fileName=${i%.*}
		FileNameNoExt=$(basename ${fileName})
		newFileName="${FileNameNoExt}_HomerAnno.txt"
		annotatePeaks.pl ${i} ${genome} > ${outDir}/${newFileName}
	done < ${bedFile}
fi