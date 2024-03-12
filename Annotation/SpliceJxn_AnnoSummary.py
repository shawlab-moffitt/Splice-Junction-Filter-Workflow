#! /bin/python3.8.6

import sys
import argparse
import os.path
import pandas as pd



parser=argparse.ArgumentParser(description='Combine summary files to generate global sample summary.')

parser.add_argument('-j','--jxns', const=1, nargs='?', required=True, default=1, metavar='',help='Junction List File')
parser.add_argument('-g','--gencodereg', const=1, nargs='?', required=False, default=1, metavar='', help='Regtools gencode annotaiton file')
parser.add_argument('-e','--ensgenereg', const=1, nargs='?', required=False, default=1, metavar='', help='Regtools ensgene annotaiton file')
parser.add_argument('-n','--ncbirefseqreg', const=1, nargs='?', required=False, default=1, metavar='', help='Regtools ncbirefseq annotaiton file')
parser.add_argument('-k','--knowngenereg', const=1, nargs='?', required=False, default=1, metavar='', help='Regtools knowngene annotaiton file')
parser.add_argument('-r','--refgenereg', const=1, nargs='?', required=False, default=1, metavar='', help='Regtools refgene annotaiton file')
parser.add_argument('-A','--appris', const=1, nargs='?', required=False, default=1, metavar='', help='Appris from BED nucleotide annotaiton file')
parser.add_argument('-a','--appris20', const=1, nargs='?', required=False, default=1, metavar='', help='Appris 20 nucleotide annotaiton file')
parser.add_argument('-b','--appris40', const=1, nargs='?', required=False, default=1, metavar='', help='Appris 40 nucleotide annotaiton file')
parser.add_argument('-x','--homer20right', const=1, nargs='?', required=False, default=1, metavar='', help='Homer 20 right annotation file')
parser.add_argument('-o','--homer40right', const=1, nargs='?', required=False, default=1, metavar='', help='Homer 40 right annotation file')
parser.add_argument('-m','--homer20left', const=1, nargs='?', required=False, default=1, metavar='', help='Homer 20 left annotation file')
parser.add_argument('-R','--homer40left', const=1, nargs='?', required=False, default=1, metavar='', help='Homer 40 left annotation file')
parser.add_argument('-O','--output', const=1, nargs='?', required=False, default=1, metavar='', help='Output file name')

args=parser.parse_args()

####----Input----####
jxnList=args.jxns
regGencodeFile=args.gencodereg
regEnsgeneFile=args.ensgenereg
regNCBIrefSeqFile=args.ncbirefseqreg
regKnowngeneFile=args.knowngenereg
regRefgeneFile=args.refgenereg
apprisFile=args.appris
appris20File=args.appris20
appris40File=args.appris40
homer20rFile=args.homer20right
homer40rFile=args.homer40right
homer20lFile=args.homer20left
homer40lFile=args.homer40left
outFile=args.output


####----Junction List File----####

jxn = pd.read_csv(jxnList, sep = '\t', header = 0)
jxn['Junction'] = jxn.iloc[:, 0]
first_column = jxn.pop('Junction')
jxn.insert(0, 'Junction', first_column)

SummaryDF = jxn


####----Regtools Files----####

if regGencodeFile == 1:
	print("Regtools gencode annotaiton file missing")
elif os.path.isfile(regGencodeFile) == False:
	print("Regtools gencode annotaiton file missing")
else: 
	regGencode = pd.read_csv(regGencodeFile, sep = '\t', header = 0)
	regGencode = regGencode.add_suffix('_Regtools_Gencode')
	regGencode = regGencode.rename({'name_Regtools_Gencode': 'Junction'}, axis=1)
	SummaryDF = pd.merge(SummaryDF,regGencode,on="Junction", how="outer")

if regEnsgeneFile == 1:
	print("Regtools ensGene annotaiton file missing")
elif os.path.isfile(regEnsgeneFile) == False:
	print("Regtools ensGene annotaiton file missing")
else: 
	regEnsgene = pd.read_csv(regEnsgeneFile, sep = '\t', header = 0)
	regEnsgene = regEnsgene.add_suffix('_Regtools_ensGene')
	regEnsgene = regEnsgene.rename({'name_Regtools_ensGene': 'Junction'}, axis=1)
	SummaryDF = pd.merge(SummaryDF,regEnsgene,on="Junction", how="outer")

if regNCBIrefSeqFile == 1:
	print("Regtools NCBIrefSeq annotaiton file missing")
elif os.path.isfile(regNCBIrefSeqFile) == False:
	print("Regtools NCBIrefSeq annotaiton file missing")
else:
	regNCBIrefSeq = pd.read_csv(regNCBIrefSeqFile, sep = '\t', header = 0)
	regNCBIrefSeq = regNCBIrefSeq.add_suffix('_Regtools_NCBIrefSeq')
	regNCBIrefSeq = regNCBIrefSeq.rename({'name_Regtools_NCBIrefSeq': 'Junction'}, axis=1)
	SummaryDF = pd.merge(SummaryDF,regNCBIrefSeq,on="Junction", how="outer")

if regKnowngeneFile == 1:
	print("Regtools Knowngene annotaiton file missing")
elif os.path.isfile(regKnowngeneFile) == False:
	print("Regtools Knowngene annotaiton file missing")
else: 
	regKnowngene = pd.read_csv(regKnowngeneFile, sep = '\t', header = 0)
	regKnowngene = regKnowngene.add_suffix('_Regtools_knownGene')
	regKnowngene = regKnowngene.rename({'name_Regtools_knownGene': 'Junction'}, axis=1)
	SummaryDF = pd.merge(SummaryDF,regKnowngene,on="Junction", how="outer")

if regRefgeneFile == 1:
	print("Regtools Refgene annotaiton file missing")
elif os.path.isfile(regRefgeneFile) == False:
	print("Regtools Refgene annotaiton file missing")
else: 
	regRefgene = pd.read_csv(regRefgeneFile, sep = '\t', header = 0)
	regRefgene = regRefgene.add_suffix('_Regtools_refGene')
	regRefgene = regRefgene.rename({'name_Regtools_refGene': 'Junction'}, axis=1)
	SummaryDF = pd.merge(SummaryDF,regRefgene,on="Junction", how="outer")


####----Appris Files----####

if apprisFile == 1:
	print("Appris annotation from BED file is missing")
elif os.path.isfile(apprisFile) == False:
	print("Appris annotation from BED file is missing")
else: 
	appris = pd.read_csv(apprisFile, sep = '\t', header = None, names = ["Chr","Start","Stop","Junction","Score","Strand","Left_Exon","Right_Exon","Skipped_Exons"])
	appris = appris.add_suffix('_Appris')
	appris = appris.rename({'Junction_Appris': 'Junction'}, axis=1)
	SummaryDF = pd.merge(SummaryDF,appris,on="Junction", how="outer")

if appris20File == 1:
	print("Regtools Appris_20 annotaiton file missing")
elif os.path.isfile(appris20File) == False:
	print("Regtools Refgene annotaiton file missing")
else: 
	appris20 = pd.read_csv(appris20File, sep = '\t', header = 0)
	appris20 = appris20.add_suffix('_Appris_20')
	appris20 = appris20.rename({'Junction_Appris_20': 'Junction'}, axis=1)
	SummaryDF = pd.merge(SummaryDF,appris20,on="Junction", how="outer")

if appris40File == 1:
	print("Regtools Appris_40 annotaiton file missing")
elif os.path.isfile(appris40File) == False:
	print("Regtools Appris_40 annotaiton file missing")
else: 
	appris40 = pd.read_csv(appris40File, sep = '\t', header = 0)
	appris40 = appris40.add_suffix('_Appris_40')
	appris40 = appris40.rename({'Junction_Appris_40': 'Junction'}, axis=1)
	SummaryDF = pd.merge(SummaryDF,appris40,on="Junction", how="outer")



####----Homer Files----####

if homer20rFile == 1:
	print("Regtools Homer_20Right annotaiton file missing")
elif os.path.isfile(homer20rFile) == False:
	print("Regtools Homer_20Right annotaiton file missing")
else: 
	homer20r = pd.read_csv(homer20rFile, sep = '\t', header = 0)
	homer20r = homer20r.add_suffix('_Homer_20Right')
	homer20r = homer20r.rename(columns={homer20r.columns[0]: 'Junction'})
	SummaryDF = pd.merge(SummaryDF,homer20r,on="Junction", how="outer")

if homer40rFile == 1:
	print("Regtools Homer_40Right annotaiton file missing")
elif os.path.isfile(homer40rFile) == False:
	print("Regtools ApHomer_40Rightpris_40 annotaiton file missing")
else: 
	homer40r = pd.read_csv(homer40rFile, sep = '\t', header = 0)
	homer40r = homer40r.add_suffix('_Homer_40Right')
	homer40r = homer40r.rename(columns={homer40r.columns[0]: 'Junction'})
	SummaryDF = pd.merge(SummaryDF,homer40r,on="Junction", how="outer")

if homer20lFile == 1:
	print("Regtools Homer_20Left annotaiton file missing")
elif os.path.isfile(homer20lFile) == False:
	print("Regtools Homer_20Left annotaiton file missing")
else: 
	homer20l = pd.read_csv(homer20lFile, sep = '\t', header = 0)
	homer20l = homer20l.add_suffix('_Homer_20Left')
	homer20l = homer20l.rename(columns={homer20l.columns[0]: 'Junction'})
	SummaryDF = pd.merge(SummaryDF,homer20l,on="Junction", how="outer")

if homer40lFile == 1:
	print("Regtools Homer_40Left annotaiton file missing")
elif os.path.isfile(homer40lFile) == False:
	print("Regtools Homer_40Left annotaiton file missing")
else: 
	homer40l = pd.read_csv(homer40lFile, sep = '\t', header = 0)
	homer40l = homer40l.add_suffix('_Homer_40Left')
	homer40l = homer40l.rename(columns={homer40l.columns[0]: 'Junction'})
	SummaryDF = pd.merge(SummaryDF,homer40l,on="Junction", how="outer")


SummaryDF.to_csv(outFile, sep = '\t', index = False, na_rep = "NaN")






