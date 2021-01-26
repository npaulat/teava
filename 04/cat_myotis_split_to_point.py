import os
import sys
import argparse
import pandas as pd
import numpy as np

### Version 2, created 27 August 2020 by Nicole S Paulat ###

### Reformats concatenated, unique-hit, headerless MELT vcf files, to return SPLIT hits to single base point insertions for accurate downstream analysis
### Also removes any hits that are less than 10kb from a scaffold end to limit any effects of misassembly on scaffold edges

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
	parser.add_argument('-in', '--index', type=str, help='Path to reference genome index file (.fai)', required=True)
	args = parser.parse_args()
	CAT = args.input
	DIR = args.directory
	OUTDIR = args.outdir
	INDEX = args.index

	return CAT, DIR, OUTDIR, INDEX

CAT, DIR, OUTDIR, INDEX = get_args()

if DIR == ".":
	DIR = os.getcwd()
if OUTDIR == ".":
	OUTDIR = os.getcwd()

BASENAME = os.path.basename(CAT).split("_dupsless")[0]
INPUT = os.path.join(DIR, CAT)
OUTBASE = os.path.join(OUTDIR, BASENAME)

HEADERS = ['#CHROM', 'POS', 'END', 'ASSESS', 'SVTYPE', 'SVLENGTH', 'ORIENTATION', 'Austroriparius', 'Brandtii', 'Ciliolabrum', 'Davidii', 'Myotis', 'Occultus', 'Sept_TTU', 'Thysanodes', 'Velifer', 'Vivesi', 'Yumanensis', 'FAMILY', 'MODULE']

DF = pd.read_csv(INPUT, sep='\t', names=HEADERS)
DF['END'] = np.where((DF['MODULE'] == 'SPLIT'), DF['POS'], DF['END'])
DF = DF[DF['POS'] >= 10000]
with open(INDEX, "r") as i:
	SCAFFOLD_DICT = {}
	for LINE in i:
		RECORD = LINE.split('\t')
		SCAFFOLD = RECORD[0]
		SCAFF_LEN = int(RECORD[1])
		try:
			SCAFFOLD_DICT[RECORD[0]].append(SCAFF_LEN)
		except KeyError:
			SCAFFOLD_DICT[SCAFFOLD] = [SCAFF_LEN]
for SCAFF_INFO in SCAFFOLD_DICT:
	DF['END'] = np.where(((DF['#CHROM'] == SCAFF_INFO) & (DF['END'] >= (SCAFFOLD_DICT[SCAFF_INFO][0] - 10000))), 0, DF['END'])
DF = DF[DF['END'] != 0]
DF = DF.sort_values(by=['#CHROM', 'POS'], ascending=[True, True])
DF.drop(columns=['ASSESS', 'FAMILY'], inplace=True)
DF.to_csv(OUTBASE+"_final_tosort.bed", sep='\t', header=True, index=False)
DF.to_csv(OUTBASE+"_split_point.bed", sep='\t', header=False, index=False)
