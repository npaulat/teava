import os
import sys
import argparse
import pandas as pd

### Version 1, created 25 August 2020 by Nicole S Paulat ###

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

BASENAME = os.path.basename(RAW_CAT).split("_unfiltered")[0]
MUT_RATE = "m" + str(MUTRATE)
INPUT = os.path.join(DIR, RAW_CAT)
OUTBASE = os.path.join(OUTDIR, BASENAME)

SPECIES_LIST = ['Austroriparius', 'Brandtii', 'Ciliolabrum', 'Davidii', 'Occultus', 'Sept_TTU', 'Thysanodes', 'Velifer', 'Vivesi', 'Yumanensis']

if MODULE == 'split':
	HEADERS1 = ['#CHROM', 'POS', 'END', 'FILTER', 'ASSESS', 'SVTYPE', 'SVLENGTH', 'ORIENTATION', 'Austroriparius', 'Brandtii', 'Ciliolabrum', 'Davidii', 'Occultus', 'Sept_TTU', 'Thysanodes', 'Velifer', 'Vivesi', 'Yumanensis', 'MODULE']
	DF = pd.read_csv(INPUT, sep='\t', names=HEADERS1)
	DF = DF[DF.FILTER == 'PASS']
	DF.to_csv(OUTBASE+"_pass.bed", sep='\t', header=True, index=False)
	DF.drop(columns=['FILTER'], inplace=True)
	DF = DF[DF['ASSESS'] >= 2]
	DF.to_csv(OUTBASE+"_assess2.bed", sep='\t', header=True, index=False)
	for SPECIES in SPECIES_LIST:
		DF = DF[DF[SPECIES] != "./."]
	DF.to_csv(OUTBASE+"_filtered_headers.bed", sep='\t', header=True, index=False)
	DF.to_csv(OUTBASE+"_filtered.bed", sep='\t', header=False, index=False)
elif MODULE == 'del':
	HEADERS2 = ['#CHROM', 'POS', 'END', 'SVTYPE', 'SVLENGTH', 'Austroriparius', 'Brandtii', 'Ciliolabrum', 'Davidii', 'Occultus', 'Sept_TTU', 'Thysanodes', 'Velifer', 'Vivesi', 'Yumanensis', 'MODULE']
	DF = pd.read_csv(INPUT, sep='\t', names=HEADERS2)
	for SPECIES in SPECIES_LIST:
		DF = DF[DF[SPECIES] != "./."]
	DF.to_csv(OUTBASE+"_filtered_headers.bed", sep='\t', header=True, index=False)
	DF.to_csv(OUTBASE+"_filtered.bed", sep='\t', header=False, index=False)
else:
	sys.exit('Invalid option; -m option must be "split" or "del"')
