#! /bin/python3.8.6

import sys
import argparse

parser=argparse.ArgumentParser(description='Annotate junctions with Appris Exon Annotation.')

parser.add_argument('-i','--input', const=1, nargs='?', required=True, default=1, metavar='',help='Junction input file')
parser.add_argument('-b','--bed', required=False,help='Flag if input file juncitons are in bed file format', action = 'store_true')
parser.add_argument('-header','--header', required=False,help='Flag if input file has header', action = 'store_true')
parser.add_argument('-a','--anno', const=1, nargs='?', required=True, default=1, metavar='',help='Annotation file')
parser.add_argument('-o','--output', const=1, nargs='?', required=False, default=1, metavar='', help='Output file name')


args=parser.parse_args()

in_file = args.input
bed_flag = args.bed
header_flag = args.header
anno_file = args.anno
out_file = args.output

outFile=open(out_file, 'w')


## Make annotation dictionary
anno_dict={}
with open(anno_file,"r") as anno:
	for line in anno:
		columns=line.strip('\n').split('\t')
		gene=columns[0].split("_")[0]
		CHR=columns[0].split("_")[1]
		start=columns[0].split("_")[2]
		stop=columns[0].split("_")[3]
		strand=columns[0].split("_")[4]
		pos='-'.join([start,stop])
		region=':'.join([CHR,pos,strand])
		anno_dict[region] = ["_".join([gene,region,columns[1]])]


with open(in_file,"r") as infile:
	if header_flag==True:
		header_line=infile.readline().strip('\n').split('\t')
		outFile.write('\t'.join(header_line+["Left_Exon","Right_Exon","Skipped_Exon"])+'\n')
	for line in infile:
		columns=line.strip('\n').split('\t')
		if bed_flag==True:
			## Break down junction data
			CHR=columns[0]
			start=columns[1]
			stop=columns[2]
			strand=columns[5]
		else:
			## Assume first column is in chr#:start-stop:strand format
			jxn=columns[0]
			CHR=jxn.split(":")[0]
			start=jxn.split(":")[1].split("-")[0]
			stop=jxn.split(":")[1].split("-")[1]
			strand=jxn.split(":")[2]
		## Get junction regions
		startReg=[(int(start)-5),(int(start)+5)]
		stopReg=[(int(stop)-5),(int(stop)+5)]
		totalReg=[startReg[0],stopReg[1]]
		## Subset annotaiton dictionary for correct chromosome and strand
		anno_val_sub=[val for key, val in anno_dict.items() if key.startswith(CHR+":") and key.endswith(":"+strand)]
		anno_val_sub=[item for sublist in anno_val_sub for item in sublist]
		anno_key_sub=[key for key, val in anno_dict.items() if key.startswith(CHR+":") and key.endswith(":"+strand)]
		## Get indices for exons that satisfy conditions
		# Upstream
		#startPos=[i.split(':')[1].split('-')[0] for i in anno_key_sub]
		#startInd=[i for i in range(len(startPos)) if startReg[0] <= int(startPos[i]) <= startReg[1]]
		stopPos=[i.split(':')[1].split('-')[1] for i in anno_key_sub]
		startInd=[i for i in range(len(stopPos)) if startReg[0] <= int(stopPos[i]) <= startReg[1]]
		# Downstream
		#stopPos=[i.split(':')[1].split('-')[1] for i in anno_key_sub]
		#stopInd=[i for i in range(len(stopPos)) if stopReg[0] <= int(stopPos[i]) <= stopReg[1]]
		startPos=[i.split(':')[1].split('-')[0] for i in anno_key_sub]
		stopInd=[i for i in range(len(startPos)) if stopReg[0] <= int(startPos[i]) <= stopReg[1]]
		# Skipped
		#startSkipInd=[i for i in range(len(startPos)) if startReg[0] >= int(startPos[i])]
		#stopSkipInd=[i for i in range(len(stopPos)) if stopReg[1] <= int(stopPos[i])]
		startSkipInd=[i for i in range(len(startPos)) if startReg[0] <= int(startPos[i])]
		stopSkipInd=[i for i in range(len(stopPos)) if stopReg[1] >= int(stopPos[i])]
		skipInd=list(set(startSkipInd) & set(stopSkipInd))
		## Get exons based off indices
		startExons=','.join([anno_val_sub[i] for i in startInd])
		if (len(startExons) == 0):
			startExons="NA"
		stopExons=','.join([anno_val_sub[i] for i in stopInd])
		if (len(stopExons) == 0):
			stopExons="NA"
		skipExons=','.join([anno_val_sub[i] for i in skipInd])
		if (len(skipExons) == 0):
			skipExons="NA"
		outFile.write('\t'.join(columns+[startExons,stopExons,skipExons])+'\n')
		

outFile.close()






