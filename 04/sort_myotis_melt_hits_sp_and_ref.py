import sys
import os
import argparse
import itertools
import numpy as np
import pandas as pd

### Version 2 created 12 August 2020 by Nicole S Paulat ###
### Arguments: input file of filtered MELT hits, output dir, mutation rate

def get_args():
	parser = argparse.ArgumentParser(description="Filter and sort differential MELT results by SpeciesA x M. myotis per element in each TE superfamily >> ouput files of sorted coordinates per TE superfamily", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	#required = parser.add_argument_group('required arguments')
	parser.add_argument('-i', '--input', type=str, help='File of filtered and reformatted MELT hits', required=True)
	parser.add_argument('-m', '--rate', type=int, help="Maximum mutation rate used to make the MELT hits file, e.g. 8, 10, or 15", required=True)
	parser.add_argument('-od', '--outdir', type=str, help='Location of directory for output files', default=".")

	args = parser.parse_args()
	MELT_HITS = args.input
	RATE = args.rate
	OUTDIR = args.outdir
	
	return MELT_HITS, RATE, OUTDIR
	
MELT_HITS, RATE, OUTDIR = get_args()
if OUTDIR == '.':
	OUTDIR = os.getcwd()

MUT_RATE = "mut" + str(RATE)

### Get full path of input file
INPUT = os.path.join(OUTDIR, MELT_HITS)

INPUT_NAME = os.path.basename(INPUT).split(".")[0]
OUT_SUBDIR = os.path.join(OUTDIR, INPUT_NAME)
if os.path.exists(OUT_SUBDIR):
	print("Input subdirectory {} exists.".format(OUT_SUBDIR))
else:
	os.mkdir(OUT_SUBDIR)
os.chdir(OUT_SUBDIR)


### Make dataframe of the input file, headers as first line by default
df = pd.read_csv(INPUT, sep='\t')

### Create summary file
SUM_FILE = INPUT_NAME + "_groups_summary_stats.txt"
SUM_OUTPUT = os.path.join(OUT_SUBDIR, SUM_FILE)
s = open(SUM_OUTPUT, "w+")
s.write("Summary of final MELT hit stats for the {} dataset\n\n".format(MUT_RATE))

###Define species of interest for comparisons of Species A vs M. Myotis (reference)
SPECIES_LIST = ['Austroriparius', 'Brandtii', 'Ciliolabrum', 'Davidii', 'Occultus', 'Sept_TTU', 'Thysanodes', 'Velifer', 'Vivesi', 'Yumanensis']

### INCLUDE TE name check - compare list of TE names in script to those in the file
### If any in the file are not within the script list, throw warning
### If any in script list missing from file, throw warning

TE_LIST = ['HAL1-1_mMyo', 'HAL1-1A_ML', 'HAL1-1B_ML', 'HAL1-1E_ML', 'L1MAB_ML', 'L1MAB2_ML', 'Ves', 'Ves1', 'Ves10', 'Ves11', 'Ves12', 'Ves13', 'Ves17', 'Ves2_mMyo', 'Ves20', 'Ves21', 'Ves25', 'Ves27', 'Ves3', 'Ves3_ML', 'Ves31', 'Ves35', 'Ves37', 'Ves38', 'Ves39', 'hAT1_ML', 'hAT1_mMyo', 'hAT2_ML', 'hat3_ML', 'nhAT_186_ML', 'nhAT17_ML', 'nhAT1a_ML', 'nhAT1b_ML', 'nhAT2_730_ML', 'nhAT2_ML', 'nhAT-3_EF', 'nhAT3_mMyo', 'nhAT34_ML', 'nhAT37_ML', 'nhAT6_ML', 'nhAT70_ML', 'SPIN_Ml', 'SPIN_NA_1_Et', 'SPIN_NA_10_Ml', 'SPIN_NA_7_Ml', 'SPIN_NA_8_Ml', 'SPIN_NA_9_Ml', 'Myotis_piggyBac', 'npiggy_156_ML', 'npiggy1_mMyo', 'npiggy111_ML', 'npiggy165_ML', 'npiggy2_345_ML', 'npiggy2_41_ML', 'npiggy2_mMyo', 'npiggy259_ML', 'npiggy269a_ML', 'npiggy3_Mlyr', 'npiggy4_mMyo', 'piggyBac_2a_Mm', 'piggyBac2_ML', 'piggyBac2_Mm', 'piggyBac2b_Mm', 'MARIN1_ML', 'Mariner1_ML', 'Mariner3_Ml', 'Mariner-84_Hsal', 'Mariner-84a_Hsal', 'Myotis_Tc1_ML', 'Myotis_Tc2_ML', 'nMar382_Ml', 'nMariner-5_EF', 'nMariner-7_EF', 'Tc1_94_ML', 'Tc2_122_ML', 'BAR1_ML', 'HeliBat_N1a_ML', 'HeliBat_N1b_ML', 'HeliBat_N1c_ML', 'HeliBat1', 'HeliBatN2', 'Helitron_R25_ML', 'Helitron10_Myo', 'Helitron-12_EF', 'Helitron12_mMyo', 'Helitron13_mMyo', 'Helitron-14_EF', 'Helitron15_mMyo', 'Helitron18_mMyo', 'Helitron-19_EF', 'Helitron19_mMyo', 'Helitron2_mMyo', 'Helitron20_mMyo', 'Helitron22_mMyo', 'Helitron-26_EF', 'Helitron26_mMyo', 'Helitron-27_EF', 'Helitron29_mMyo', 'Helitron3_mMyo', 'Helitron-30_EF', 'Helitron32_mMyo', 'Helitron-33_EF', 'Helitron-38_EF', 'Helitron-4_EF', 'Helitron-43_EF', 'Helitron5_mMyo', 'Helitron-63_EF', 'Helitron7_mMyo', 'Helitron8_mMyo', 'Helitron-9_EF']

TYPE_DATA = df['SVTYPE'].unique()
TYPE_LIST = np.array(TYPE_DATA).tolist()
TYPE_LIST.sort(key=lambda x: x.lower())

for SV_TYPE in TYPE_LIST:
	if SV_TYPE not in TE_LIST:
		sys.exit("{} is not in TE sublibrary list!!".format(SV_TYPE))

### Define TE superfamilies and their associated (relevant) TE subfamilies as a dictionary
### Current TE superfamilies = HAL, hAT, Helitron, L1, piggyBac, TcMariner, Ves (SINE)

TE_DICT = {
	"HAL": ["HAL1-1_mMyo", "HAL1-1A_ML", "HAL1-1B_ML", "HAL1-1E_ML"],
	"LINE": ["L1MAB_ML", "L1MAB2_ML"],
	"SINE": ["Ves", "Ves1", "Ves10", "Ves11", "Ves12", "Ves13", "Ves17", "Ves2_mMyo", "Ves20", "Ves21", "Ves25", "Ves27", "Ves3", "Ves3_ML", "Ves31", "Ves35", "Ves37", "Ves38", "Ves39"],
	"hAT": ["hAT1_ML", "hAT1_mMyo", "hAT2_ML", "hat3_ML", "nhAT_186_ML", "nhAT17_ML", "nhAT1a_ML", "nhAT1b_ML", "nhAT2_730_ML", "nhAT2_ML", "nhAT-3_EF", "nhAT3_mMyo", "nhAT34_ML", "nhAT37_ML", "nhAT6_ML", "nhAT70_ML", "SPIN_Ml", "SPIN_NA_1_Et", "SPIN_NA_10_Ml", "SPIN_NA_7_Ml", "SPIN_NA_8_Ml", "SPIN_NA_9_Ml"],
	"PiggyBac": ["Myotis_piggyBac", "npiggy_156_ML", "npiggy1_mMyo", "npiggy111_ML", "npiggy165_ML", "npiggy2_345_ML", "npiggy2_41_ML", "npiggy2_mMyo", "npiggy259_ML", "npiggy269a_ML", "npiggy3_Mlyr", "npiggy4_mMyo", "piggyBac_2a_Mm", "piggyBac2_ML", "piggyBac2_Mm", "piggyBac2b_Mm"],
	"TcMariner": ["MARIN1_ML", "Mariner1_ML", "Mariner3_Ml", "Mariner-84_Hsal", "Mariner-84a_Hsal", "Myotis_Tc1_ML", "Myotis_Tc2_ML", "nMar382_Ml", "nMariner-5_EF", "nMariner-7_EF", "Tc1_94_ML", "Tc2_122_ML"],
	"Helitron": ["BAR1_ML", "HeliBat_N1a_ML", "HeliBat_N1b_ML", "HeliBat_N1c_ML", "HeliBat1", "HeliBatN2", "Helitron_R25_ML", "Helitron10_Myo", "Helitron-12_EF", "Helitron12_mMyo", "Helitron13_mMyo", "Helitron-14_EF", "Helitron15_mMyo", "Helitron18_mMyo", "Helitron-19_EF", "Helitron19_mMyo", "Helitron2_mMyo", "Helitron20_mMyo", "Helitron22_mMyo", "Helitron-26_EF", "Helitron26_mMyo", "Helitron-27_EF", "Helitron29_mMyo", "Helitron3_mMyo", "Helitron-30_EF", "Helitron32_mMyo", "Helitron-33_EF", "Helitron-38_EF", "Helitron-4_EF", "Helitron-43_EF", "Helitron5_mMyo", "Helitron-63_EF", "Helitron7_mMyo", "Helitron8_mMyo", "Helitron-9_EF"]
};

### Extract differential hits that match each species pair and each TE superfamily
### Differential hit = Species A : 0/1 or 1/1 vs M. Myotis : 0/0; or vice versa
### Write all species-pair hits for all elements w/n each superfamily to appropriate output file

# For each species from the list:
for SPECIES in SPECIES_LIST:
	s.write("{} TE hit summary:\n".format(SPECIES))
	# For each TE superfamily in the TE dictionary:
	for TE_FAMILY in TE_DICT:
		s.write("{} superfamily\n".format(TE_FAMILY))
		# Make the appropriate output file given the SPECIES, TE superfamily, and mutation rate
		OUT_FILE1 = SPECIES + "_" + TE_FAMILY + "_SPLIT_" + MUT_RATE + ".bed"
		OUTPUT1 = os.path.join(OUT_SUBDIR, OUT_FILE1)
		OUT_FILE2 = SPECIES + "_" + TE_FAMILY + "_DEL_" + MUT_RATE + ".bed"
		OUTPUT2 = os.path.join(OUT_SUBDIR, OUT_FILE2)
		OUT_FILE3 = SPECIES + "_" + TE_FAMILY + "_all_" + MUT_RATE + ".bed"
		OUTPUT3 = os.path.join(OUT_SUBDIR, OUT_FILE3)
		# Initialize the output list for this species-pair and TE superfamily
		FAMILY_SP_PAIR_HITS = []
		#Initialize the SPLIT and DELETION output lists for this species-pair and TE superfamily
		FAMILY_SP_PAIR_SPLIT_HITS = []
		FAMILY_SP_PAIR_DEL_HITS = []
		# Initialize superfamily insertion count
		FAMILY_HIT_COUNT = 0
		# For each TE subfamily within the TE superfamily:
		for SUBFAMILY in TE_DICT[TE_FAMILY]:
			# Initialize the list for this species-pair and TE subfamily
			TE_SPECIES_PAIR_HITS = []
			TE_SPECIES_PAIR_SPLIT_HITS = []
			TE_SPECIES_PAIR_DEL_HITS = []
			# Make dataframes for differential hits for the species pair
			# These will correspond to subsets of SPLIT and DELETION hits, respectively
			df1 = df[(df['SVTYPE'] == SUBFAMILY) & (df[SPECIES] != "0/0") & (df['Myotis'] == "0/0")]
			df2 = df[(df['SVTYPE'] == SUBFAMILY) & (df[SPECIES] == "0/0") & (df['Myotis'] != "0/0")]
			### To simply have a single dataset of all TE insertions, use the following commands
			# Concatenate the two dataframes
			df3 = pd.concat([df1, df2])
			# Subset the differential hits dataframe to have just the coordinates
			df4 = df3[['#CHROM', 'POS', 'END']]
			# Convert the dataframe of coordinates into a list of lists
			TE_SPECIES_PAIR_HITS = df4.values.tolist()
			# Subset pf coordinates for non-reference TE insertions (SPLIT hits)
			df_SPLIT = df1[['#CHROM', 'POS', 'END']]
			df_SPLIT.reset_index(inplace=True, drop=True)
			SPLIT_COUNT = len(df_SPLIT.index)
			s.write("{} nonreference (SPLIT) insertions = {}\n".format(SUBFAMILY, str(SPLIT_COUNT)))
			# Subset of coordinates for reference TE insertions (DELETION hits)
			df_DEL = df2[['#CHROM', 'POS', 'END']]
			df_DEL.reset_index(inplace=True, drop=True)
			DEL_COUNT = len(df_DEL.index)
			s.write("{} reference (DELETION) insertions = {}\n".format(SUBFAMILY, str(DEL_COUNT)))
			# Convert the dataframe of coordinates into a list of lists
			TE_SPECIES_PAIR_SPLIT_HITS = df_SPLIT.values.tolist()
			# Convert the dataframe of coordinates into a list of lists
			TE_SPECIES_PAIR_DEL_HITS = df_DEL.values.tolist()
			FAMILY_HIT_COUNT = FAMILY_HIT_COUNT + (SPLIT_COUNT + DEL_COUNT)
			# For each of the hits in the TE subfamily hit list generated:
			# Append to the TE superfamily hits coordinate list
			for SUBFAM_HIT in TE_SPECIES_PAIR_HITS:
				FAMILY_SP_PAIR_HITS.append(SUBFAM_HIT)
			for SUBFAM_HIT in TE_SPECIES_PAIR_SPLIT_HITS:
				FAMILY_SP_PAIR_SPLIT_HITS.append(SUBFAM_HIT)
			for SUBFAM_HIT in TE_SPECIES_PAIR_DEL_HITS:
				FAMILY_SP_PAIR_DEL_HITS.append(SUBFAM_HIT)
		# Once all TE subfamilies for a given superfamily are appended, sort by CHROM and POS
		FAMILY_SP_PAIR_HITS.sort(key=lambda x:(x[0], x[1]))
		FAMILY_SP_PAIR_SPLIT_HITS.sort(key=lambda x:(x[0], x[1]))
		FAMILY_SP_PAIR_DEL_HITS.sort(key=lambda x:(x[0], x[1]))
		# Write out number of MELT hits for the TE superfamily
		s.write("{} superfamily total insertions = {}\n\n".format(TE_FAMILY, str(FAMILY_HIT_COUNT)))
		# Write out each hit within the TE superfamily to the output file
		with open(OUTPUT1, "w+") as f:
			for HIT in FAMILY_SP_PAIR_SPLIT_HITS:
				f.write("{}\t{}\t{}\n".format(HIT[0], HIT[1], HIT[2]))
		with open(OUTPUT2, "w+") as g:
			for HIT in FAMILY_SP_PAIR_DEL_HITS:
				g.write("{}\t{}\t{}\n".format(HIT[0], HIT[1], HIT[2]))
		with open(OUTPUT3, "w+") as j:
			for HIT in FAMILY_SP_PAIR_HITS:
				j.write("{}\t{}\t{}\n".format(HIT[0], HIT[1], HIT[2]))
