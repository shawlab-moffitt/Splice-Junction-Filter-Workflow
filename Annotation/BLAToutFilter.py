#! /bin/python3.8.6

import sys

inFileName=sys.argv[1]     #Samtools faidx Summary File
pslFileName=sys.argv[2]	   #BLAT PSL
blat8FileName=sys.argv[3]  #BLAT8
outFileName=sys.argv[4]    #Output file name

outFile=open(outFileName, 'w')

pslPass=[]
pslHardPass=[]
inSize=[]

with open(pslFileName,"r") as psl:
	pslLineCount=0
	for line in psl:
		if line.startswith("-"):
			for line in psl:
				columns=line.strip('\n').split('\t')
				jxn=columns[9]
				CHR=jxn.split(":")[0]
				TgapCount=columns[6]
				name=columns[13]
				size=str(int(columns[10])/2).strip(".0")
				inSize.append(columns[10])
				blockCount=columns[17]
				blockSize=columns[18]
				blockSize1=blockSize.split(",")[0]
				blockSize2=blockSize.split(",")[1]
				qStarts=columns[19]
				if (int(TgapCount)==1) and (int(blockCount)==2) and (str(blockSize)==size+","+size+",") and (str(qStarts)=="0,"+size+",") and (name == CHR):
					pslHardPass.append(jxn)
				#else:
				#	continue
				if len(blockSize.split(",")) == 3:
					if (int(TgapCount)==1) and (int(blockCount)==2) and ((int(blockSize1)+int(blockSize2))==int(columns[10])) and (str(qStarts)==("0,"+blockSize1+",")) and (name == CHR):
						pslPass.append(jxn)
					#else:
					#	continue
				#else:
				#	continue
psl.close()

blast8Pass=[]
with open(blat8FileName,"r") as blast8:
	blast8LineCount=0
	for line in blast8:
		columns=line.strip('\n').split('\t')
		jxn=columns[0]
		CHR=jxn.split(":")[0]
		name=columns[1]
		size=inSize[0]
		pident=columns[2]
		qstart=columns[6]
		qend=columns[7]
		if (pident=="100.00") and ((qstart=="1") and (qend==size)) and (name == CHR):
			blast8Pass.append(jxn)
		else:
			continue
blast8.close()

with open(inFileName,"r") as inFile:
	inFileLineCount=0
	for line in inFile:
		columns = line.strip('\n').split('\t')
		jxn=columns[0]
		if columns[0] == "Junction":
			outFile.write('\t'.join(columns+["PSL_Jxn_Relax_Filter","PSL_Jxn_Strict_Filter","Blast8_Jxn_Filter"])+'\n')
		else:
			outFile.write('\t'.join(columns+[str(pslPass.count(jxn)),str(pslHardPass.count(jxn)),str(blast8Pass.count(jxn))])+'\n')

inFile.close()
outFile.close()

				