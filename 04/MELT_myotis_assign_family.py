import os
import sys
import argparse
import pandas as pd
import numpy as np

### Version 1, created 1 September 2020 by Nicole S Paulat ###


def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="Assignment of TE family/superfamily to final MELT TE BED file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	required = parser.add_argument_group('required arguments')
	#Give input coordinate file of formatted MELT SPLIT hits
	parser.add_argument('-i', '--input', help='filtered, final MELT file', required=True)
	#Argument of directory containing the list of formatted MELT hits (need full path)
	parser.add_argument('-d', '--directory', type=str, help='Path to the directory of the input file', default=".")
	#Argument of the output directory (need full path)
	parser.add_argument('-od', '--outdir', type=str, help='Location of directory for the output file', default=".")
	#Argument of the output directory (need full path)
	parser.add_argument('-tecat', '--tecategories', type=str, help='Path to the list of RepeatMasker TE names, with TE family category (te_lib_categories.txt)', required=True)
	
	args = parser.parse_args()
	INPUT = args.input
	DIR = args.directory
	OUTDIR = args.outdir
	CAT_TE_NAMES = args.tecategories

	return INPUT, DIR, OUTDIR, CAT_TE_NAMES

INPUT, DIR, OUTDIR, CAT_TE_NAMES = get_args()

if DIR == ".":
	DIR = os.getcwd()
if OUTDIR == ".":
	OUTDIR = os.getcwd()

BASENAME = os.path.basename(INPUT).split("_cat")[0]
CAT_HITS = os.path.join(DIR, INPUT)
OUTBASE = os.path.join(OUTDIR, BASENAME)
OUTPUT = OUTBASE + "_cat_final_type.bed"

HEADERS = ['#CHROM', 'POS', 'END', 'SVTYPE', 'SVLENGTH', 'ORIENTATION', 'Austroriparius', 'Brandtii', 'Ciliolabrum', 'Davidii', 'Myotis', 'Occultus', 'Sept_TTU', 'Thysanodes', 'Velifer', 'Vivesi', 'Yumanensis', 'MODULE']

CAT_DF = pd.read_csv(CAT_HITS, sep='\t')
CAT_DF.sort_values(by=['#CHROM', 'POS'], ascending=[True, True], inplace=True)
CAT_DF.insert(18, 'TYPE', '')
for LINE in open(CAT_TE_NAMES, "r"):
	CAT_NAME, GROUP = LINE.rstrip().split(" ", 1)
	TES = GROUP.split(" ")
	for TE in TES:
		CAT_DF['TYPE'] = np.where((CAT_DF['SVTYPE'] == TE), CAT_NAME, CAT_DF['TYPE'])

CAT_DF.to_csv(OUTPUT, sep='\t', header=True, index=False)
