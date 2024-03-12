
#!/bin/bash

bedFiles=$1
faFile=$2
gtfFile=$3
annoName=$4
outDir=$5

ml RegTools/0.5.2-foss-2021b

mkdir -p ${outDir}

bedCols1=$(awk -F'\t' '{print NF; exit}' ${bedFiles})
base=1
if [ "${bedCols1}" -gt "${base}" ]
then
	bedCols=$(awk -F'\t' '{print NF; exit}' ${bedFiles})
	FullFileName=${bedFiles%.*}
	FileNameNoExt=$(basename ${FullFileName})
	bed6=6
	bed12=12
	if [ "${bedCols}" -eq "${bed6}" ]
	then
		bed62bed12.pl -i ${bedFiles} > ${bedFiles%.*}.bed12
		cat ${bedFiles%.*}.bed12 | awk -F'\t' '!/_/ && !/^chrEBV/ { print }' > ${bedFiles%.*}.knownPos.bed12
		regtools junctions annotate -o ${outDir}/${FileNameNoExt}.knownPos.${annoName}_Anno.bed12 ${bedFiles%.*}.knownPos.bed12 ${faFile} ${gtfFile}
	elif [ "${bedCols}" -eq "${bed12}" ]
	then
		cat ${bedFiles%.*}.bed12 | awk -F'\t' '!/_/ && !/^chrEBV/ { print }' > ${bedFiles%.*}.knownPos.bed12
		regtools junctions annotate -o ${outDir}/${FileNameNoExt}.knownPos.${annoName}_Anno.bed12 ${bedFiles%.*}.knownPos.bed12 ${faFile} ${gtfFile}
	else
		echo "incorrect number of columns in bed file"
	fi
elif [ "${bedCols1}" -eq "${base}" ]
then
	while read i
	do
		fileName=${i%.*}
		bedCols=$(awk -F'\t' '{print NF; exit}' ${i})
		FullFileName=${fileName%.*}
		FileNameNoExt=$(basename ${FullFileName})
		bed6=6
		bed12=12
		if [ "${bedCols}" -eq "${bed6}" ]
		then
			bed62bed12.pl -i ${i} > ${i%.*}.bed12
			cat ${i%.*}.bed12 | awk -F'\t' '!/_/ && !/^chrEBV/ { print }' > ${i%.*}.knownPos.bed12
			regtools junctions annotate -o ${outDir}/${FileNameNoExt}.knownPos.${annoName}_Anno.bed12 ${i%.*}.knownPos.bed12 ${faFile} ${gtfFile}
		elif [ "${bedCols}" -eq "${bed12}" ]
		then
			cat ${i%.*}.bed12 | awk -F'\t' '!/_/ && !/^chrEBV/ { print }' > ${i%.*}.knownPos.bed12
			regtools junctions annotate -o ${outDir}/${FileNameNoExt}.knownPos.${annoName}_Anno.bed12 ${i%.*}.knownPos.bed12 ${faFile} ${gtfFile}
		else
			echo "incorrect number of columns in bed file"
		fi
	done < ${bedFiles}
fi





