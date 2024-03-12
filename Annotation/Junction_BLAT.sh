#!bin/bash

ml SAMtools/1.16.1-GCC-11.3.0

JxnList=${1}
faFile=${2}
outFile=${3}
bp=${4}

echo -e "Junction\tStartPos_${bp}_Left\tStartPos_${bp}_Right\tStopPos_${bp}_Left\tStopPos_${bp}_Right\tOuterJxn_Region" >> ${outFile}

outfa=${outFile%.*}

while read i
do
	junction=${i}
	IFS=: read -r -a junction_array <<< "${junction}"
	chr=${junction_array[0]}
	pos=${junction_array[1]}
	IFS=- read -r -a pos_array <<< "${pos}"
	start=${pos_array[0]}
	stop=${pos_array[1]}

	startDown20=$(($((start-1)) - ${bp}))
	startUp20=$((start + ${bp}))
	stopDown20=$((stop - ${bp}))
	stopUp20=$(($((stop+1)) + ${bp}))

	startDown20region=${chr}:${startDown20}-${start}
	startUp20region=${chr}:${start}-${startUp20}
	stopDown20region=${chr}:${stopDown20}-${stop}
	stopUp20region=${chr}:${stop}-${stopUp20}

	startDown20blat=$(samtools faidx ${faFile} ${startDown20region} | sed -n '2p')
	startUp20blat=$(samtools faidx ${faFile} ${startUp20region} | sed -n '2p')
	stopDown20blat=$(samtools faidx ${faFile} ${stopDown20region} | sed -n '2p')
	stopUp20blat=$(samtools faidx ${faFile} ${stopUp20region} | sed -n '2p')

	echo ">${junction}" >> ${outfa}.fa
	echo "${startDown20blat}${stopUp20blat}" >> ${outfa}.fa

	echo -e "${junction}\t${startDown20blat}\t${startUp20blat}\t${stopDown20blat}\t${stopUp20blat}\t${startDown20blat}${stopUp20blat}" >> ${outFile}
done < ${JxnList}

