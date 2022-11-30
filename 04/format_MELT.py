import os
import sys
import argparse
import pandas as pd

### Version 1, created 14 May 2019 by Nicole S Paulat ###

### Reformats concatenated, headerless MELT vcf files, into the relevant information columns, with extraneous information/columns removed, ready to use in the duplicate-removal scripts

def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="General removal of most MELT duplicate calls and overlapping calls", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	required = parser.add_argument_group('required arguments')
	#Give input coordinate file of MELT hits (CHROM# + POS)
	parser.add_argument('-i', '--input', help='raw concatenated MELT output file', required=True)
	#Argument of directory containing the list of formatted MELT hits (need full path)
	parser.add_argument('-d', '--directory', type=str, help='Path to the directory of the input file', default=".")
	#Argument of the output directory (need full path)
	parser.add_argument('-od', '--outdir', type=str, help='Location of directory for the output file', default=".")
	#Indicate if MELT SPLIT file or MELT DELETION file
	parser.add_argument('-m', '--module', type=str, help='Indicate if MELT SPLIT file (split) or MELT DELETION file (del)', required=True)
	#Indicate if MELT SPLIT file or MELT DELETION file
	parser.add_argument('-r', '--mutrate', type=int, help='Indicate maximum mutation rate of the dataset (e.g. 8)', required=True)
	args = parser.parse_args()
	RAW_CAT = args.input
	DIR = args.directory
	OUTDIR = args.outdir
	MODULE = args.module
	MUTRATE = args.mutrate

	return RAW_CAT, DIR, OUTDIR, MODULE, MUTRATE

RAW_CAT, DIR, OUTDIR, MODULE, MUTRATE = get_args()

if DIR == ".":
	DIR = os.getcwd()
if OUTDIR == ".":
	OUTDIR = os.getcwd()

BASENAME = os.path.basename(RAW_CAT).split("_results")[0]
MUT_RATE = "m" + str(MUTRATE)
INPUT = os.path.join(DIR, RAW_CAT)
OUTBASE = os.path.join(OUTDIR, BASENAME)

HEADERS = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'Austroriparius', 'Brandtii', 'Ciliolabrum', 'Davidii', 'Occultus', 'Sept_TTU', 'Sept_USDA', 'Thysanodes', 'Velifer', 'Vivesi', 'Yumanensis']
SPECIES_LIST = ['Austroriparius', 'Brandtii', 'Ciliolabrum', 'Davidii', 'Occultus', 'Sept_TTU', 'Sept_USDA', 'Thysanodes', 'Velifer', 'Vivesi', 'Yumanensis']

DF = pd.read_csv(INPUT, sep='\t', names=HEADERS)

INFORMATION = DF['INFO'].str.split(";", expand=True)

if MODULE == 'split':
	DF.drop(columns=['ID', 'REF', 'ALT', 'QUAL', 'FORMAT'], inplace=True)
	DF = DF[DF.FILTER == 'PASS']
	DF.to_csv(OUTBASE+"_pass.bed", sep='\t', header=True, index=False)
	DF.drop(columns=['FILTER'], inplace=True)
	DF['ASSESS']=INFORMATION[1].str.split("=").str[1].apply(pd.to_numeric)
	DF['SVTYPE']=INFORMATION[3].str.split("=").str[1].str.split(MUT_RATE).str[0]
	DF['SVLENGTH']=INFORMATION[4].str.split("=").str[1].apply(pd.to_numeric)
	DF['ORIENTATION']=INFORMATION[5].str.split(",").str[-1]
	DF.drop(columns=['INFO'], inplace=True)
	#DF[['ASSESS', 'SVLENGTH']] = DF[['ASSESS', 'SVLENGTH']].apply(pd.to_numeric)
	DF['END'] = DF['POS'] + DF['SVLENGTH']
	DF['MODULE'] = 'SPLIT'
	for SPECIES in SPECIES_LIST:
		DF[SPECIES] = DF[SPECIES].str.split(":").str[0]
	DF = DF[['#CHROM', 'POS', 'END', 'ASSESS', 'SVTYPE', 'SVLENGTH', 'ORIENTATION', 'Austroriparius', 'Brandtii', 'Ciliolabrum', 'Davidii', 'Occultus', 'Sept_TTU', 'Sept_USDA', 'Thysanodes', 'Velifer', 'Vivesi', 'Yumanensis', 'MODULE']]
	DF = DF[DF['ASSESS'] >= 2]
	DF.to_csv(OUTBASE+"_assess2.bed", sep='\t', header=True, index=False)
	#DF.drop(columns=['ASSESS'], inplace=True) # gives chained warning "SettingWithCopyWarning"
	#DF = DF.drop(columns=['ASSESS'])
	for SPECIES in SPECIES_LIST:
		DF = DF[DF[SPECIES] != "./."]
	DF.insert(11, 'Lucifugus', '0/0')
	DF.to_csv(OUTBASE+"_filtered.bed", sep='\t', header=True, index=False)
	DF.to_csv(OUTBASE+"_filtered_dups.bed", sep='\t', header=False, index=False)
elif MODULE == 'del':
	DF.drop(columns=['ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'FORMAT'], inplace=True)
	INFORMATION = DF['INFO'].str.split(";", expand=True)
	DF['SVTYPE']=INFORMATION[0].str.split("=").str[1].str.split(MUT_RATE).str[0]
	DF['END']=INFORMATION[1].str.split('=').str[1].apply(pd.to_numeric)
	DF['SVLENGTH']=INFORMATION[2].str.split("=").str[1].apply(pd.to_numeric)
	DF['MODULE'] = 'DELETION'
	DF.drop(columns=['INFO'], inplace=True)
	for SPECIES in SPECIES_LIST:
		DF[SPECIES] = DF[SPECIES].str.split(":").str[0]
	DF = DF[['#CHROM', 'POS', 'END', 'SVTYPE', 'SVLENGTH', 'Austroriparius', 'Brandtii', 'Ciliolabrum', 'Davidii', 'Occultus', 'Sept_TTU', 'Sept_USDA', 'Thysanodes', 'Velifer', 'Vivesi', 'Yumanensis', 'MODULE']]
	for SPECIES in SPECIES_LIST:
		DF = DF[DF[SPECIES] != "./."]
	DF.insert(11, 'Lucifugus', '0/0')
	DF.replace({'0/0': '1/1', '1/1': '0/0'}, inplace=True)
	DF.to_csv(OUTBASE+"_filtered.bed", sep='\t', header=True, index=False)
	DF.to_csv(OUTBASE+"_filtered_dups.bed", sep='\t', header=False, index=False)
	DF2 = DF[['#CHROM', 'POS', 'END', 'SVTYPE']]
	DF2.to_csv(OUTBASE+"_add_orient.bed", sep='\t', header=False, index=False)
else:
	sys.exit('Invalid option; -m option must be "split" or "del"')
