import os
import sys
import argparse
import pandas as pd
import numpy as np

### Version 3, created 12 August 2020 by Nicole S Paulat ###

### Reformats concatenated, headerless MELT vcf files, into the relevant information columns, with extraneous information/columns removed, ready to use in the duplicate-removal scripts
### This includes renaming MELT SPLIT hits to match the original TE names from RepeatMasker and the original TE library

def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="General removal of most MELT duplicate calls and overlapping calls", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	required = parser.add_argument_group('required arguments')
	#Give input coordinate file of formatted MELT SPLIT hits
	parser.add_argument('-s', '--split', help='filtered, unique, concatenated MELT SPLIT file', required=True)
	#Argument of directory containing the formatted MELT DELETION hits
	parser.add_argument('-del', '--deletion', help='filtered, unique, concatenated MELT DELETION file', required=True)
	#Argument of directory containing the list of formatted MELT hits (need full path)
	parser.add_argument('-d', '--directory', type=str, help='Path to the directory of the input file', default=".")
	#Argument of the output directory (need full path)
	parser.add_argument('-od', '--outdir', type=str, help='Location of directory for the output file', default=".")
	#Argument of the output directory (need full path)
	parser.add_argument('-mt', '--melttes', type=str, help='Path to list of MELT ZIP basenames (MELT compatible TE names = zip_te_names.txt)', required=True)
	#Argument of the output directory (need full path)
	parser.add_argument('-lib', '--telibrary', type=str, help='Path to the list of RepeatMasker TE names (te_list.txt)', required=True)
	#Argument of the output directory (need full path)
	parser.add_argument('-tecat', '--tecategories', type=str, help='Path to the list of RepeatMasker TE names, with TE family category (te_lib_categories.txt)', required=True)
	
	args = parser.parse_args()
	SPLIT = args.split
	DELETION = args.deletion
	DIR = args.directory
	OUTDIR = args.outdir
	SPLIT_TE_NAMES = args.melttes
	TE_NAMES = args.telibrary
	CAT_TE_NAMES = args.tecategories

	return SPLIT, DELETION, DIR, OUTDIR, SPLIT_TE_NAMES, TE_NAMES, CAT_TE_NAMES

SPLIT, DELETION, DIR, OUTDIR, SPLIT_TE_NAMES, TE_NAMES, CAT_TE_NAMES = get_args()

if DIR == ".":
	DIR = os.getcwd()
if OUTDIR == ".":
	OUTDIR = os.getcwd()

BASENAME = os.path.basename(SPLIT).split("_SPLIT")[0]
SPLIT_HITS = os.path.join(DIR, SPLIT)
DEL_HITS = os.path.join(DIR, DELETION)
OUTBASE = os.path.join(OUTDIR, BASENAME)
OUTPUT1 = OUTBASE + "_cat_assess_dups_headers.bed"
OUTPUT2 = OUTBASE + "_cat_assess_dups.bed"

HEADERS1 = ['#CHROM', 'POS', 'END', 'ASSESS', 'SVTYPE', 'SVLENGTH', 'ORIENTATION', 'Austroriparius', 'Brandtii', 'Ciliolabrum', 'Davidii', 'Myotis', 'Occultus', 'Sept_TTU', 'Thysanodes', 'Velifer', 'Vivesi', 'Yumanensis', 'MODULE']
HEADERS2 = ['#CHROM', 'POS', 'END', 'SVTYPE', 'SVLENGTH', 'ORIENTATION', 'Austroriparius', 'Brandtii', 'Ciliolabrum', 'Davidii', 'Myotis', 'Occultus', 'Sept_TTU', 'Thysanodes', 'Velifer', 'Vivesi', 'Yumanensis', 'MODULE']

SPLIT_DF = pd.read_csv(SPLIT_HITS, sep='\t', names=HEADERS1)
DEL_DF = pd.read_csv(DEL_HITS, sep='\t', names=HEADERS2)

DEL_DF.insert(3,'ASSESS', '4')
#SPLIT_DF.drop(columns=['ASSESS'], inplace=True)
CAT_DF = pd.concat([SPLIT_DF, DEL_DF])
CAT_DF.sort_values(by=['#CHROM', 'POS'], ascending=[True, True], inplace=True)
with open(SPLIT_TE_NAMES, "r") as d:
	SPLIT_TES = d.read().split(" ")
with open(TE_NAMES, "r") as d:
	REF_TES = d.read().split(" ")
for SPLIT_TE, REF_TE in zip(SPLIT_TES, REF_TES):
	CAT_DF['SVTYPE'] = np.where((CAT_DF['SVTYPE'] == SPLIT_TE), REF_TE, CAT_DF['SVTYPE'])
CAT_DF.insert(18, 'TYPE', '')
for LINE in open(CAT_TE_NAMES, "r"):
	CAT_NAME, GROUP = LINE.rstrip().split(" ", 1)
	TES = GROUP.split(" ")
	for TE in TES:
		CAT_DF['TYPE'] = np.where((CAT_DF['SVTYPE'] == TE), CAT_NAME, CAT_DF['TYPE'])

CAT_DF.to_csv(OUTPUT1, sep='\t', header=True, index=False)
#CAT_DF.to_csv(OUTPUT1, sep='\t', header=False, index=False)
CAT_DF.to_csv(OUTPUT2, sep='\t', header=False, index=False)
