import sys
import os
import argparse
import itertools
import numpy as np
import pandas as pd

#Arguments: input file of filtered MELT hits, output dir, mutation rate

def get_args():
	parser = argparse.ArgumentParser(description="Filter and sort differential MELT results by SpeciesA x Lucifugus per element in each TE superfamily >> ouput files of sorted coordinates per TE superfamily", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	#required = parser.add_argument_group('required arguments')
	parser.add_argument('-i', '--input', type=str, help='File of filtered and reformatted MELT hits', required=True)
	parser.add_argument('-m', '--rate', type=int, help="Maximum mutation rate used to make the MELT hits file, e.g. 5, 8, or 15", required=True)
	parser.add_argument('-od', '--outdir', type=str, help='Location of directory for output files', required=True)

	args = parser.parse_args()
	MELT_HITS = args.input
	RATE = args.rate
	OUTDIR = args.outdir
	
	return MELT_HITS, RATE, OUTDIR
	
MELT_HITS, RATE, OUTDIR = get_args()

MUT_RATE = "mut" + str(RATE)

### Get full path of input file

INPUT = os.path.join(OUTDIR, MELT_HITS)

### Make dataframe of the input file, headers as first line by default

df = pd.read_csv(INPUT, sep='\t')

###Define species of interest for comparisons of Species A vs Lucifugus (reference)

SPECIES_LIST = ['Austroriparius', 'Brandtii', 'Ciliolabrum', 'Davidii', 'Occultus', 'Sept_TTU', 'Sept_USDA', 'Thysanodes', 'Velifer', 'Vivesi', 'Yumanensis']

### INCLUDE TE name check - compare list of TE names in script to those in the file
### If any in the file are not within the script list, throw warning
### If any in script list missing from file, throw warning

TE_LIST = ['HAL1-1B_ML', 'HAL1-1E_ML', 'hAT-2N1_ML', 'hAT-2N2_ML', 'HeliBat_N1a_ML', 'HeliBat_N1b_ML', 'HeliBat_N1c_ML', 'Helitron-27_EF', 'Helitron-33_EF', 'Helitron-36_EF', 'Helitron-37_EF', 'Helitron-39B_EF', 'Helitron-3_EF', 'Helitron-52_EF', 'Helitron-54_EF', 'Helitron-63_EF', 'Helitron-8_EF', 'Helitron-9_EF', 'Helitron_R25_ML', 'L1MAB2_ML', 'L1MAB_ML', 'MARIN1_ML', 'Mariner1_ML', 'Mariner3_Ml', 'nhAT-1_EF', 'nhAT114_ML', 'nhAT171_ML', 'nhAT17_ML', 'nhAT1_ML', 'nhAT2_730_ML', 'nhAT34_ML', 'nhAT37_ML', 'nhAT3_ML', 'nhAT68_ML', 'nhAT6_ML', 'nhAT70_ML', 'nhAT_186_ML', 'nMar.Tc1_311_ML', 'nMar382_Ml', 'nMar91_ML', 'nMariner-1_EF', 'npiggy165_ML', 'npiggy259_ML', 'npiggy2_41_ML', 'npiggy_156_ML', 'nTc1_135_ML', 'nTc1_27_ML', 'piggyBac2_ML', 'SPIN_Ml', 'Tc1_94_ML', 'Ves1', 'Ves12', 'Ves13', 'Ves17', 'Ves20', 'Ves21', 'Ves25', 'Ves27', 'Ves3', 'Ves31', 'Ves35', 'Ves38', 'Ves3_ML']

TYPE_DATA = df['SVTYPE'].unique()
TYPE_LIST = np.array(TYPE_DATA).tolist()
TYPE_LIST.sort(key=lambda x: x.lower())

for SV_TYPE in TYPE_LIST:
	if SV_TYPE not in TE_LIST:
		sys.exit("{} is not in TE sublibrary list!!".format(SV_TYPE))

### Define TE superfamilies and their associated (relevant) TE subfamilies as a dictionary
### Current TE superfamilies = HAL, hAT, Helitron, L1, piggyBac, TcMariner, Ves (SINE)

TE_DICT = {
	"hal": ["HAL1-1B_ML", "HAL1-1E_ML"],
	"hat": ["hAT-2N1_ML", "hAT-2N2_ML", "nhAT_186_ML", "nhAT-1_EF", "nhAT1_ML", "nhAT114_ML", "nhAT17_ML", "nhAT171_ML", "nhAT2_730_ML", "nhAT3_ML", "nhAT34_ML", "nhAT37_ML", "nhAT6_ML", "nhAT68_ML", "nhAT70_ML", "SPIN_Ml"],
	"helitron": ["HeliBat_N1a_ML", "HeliBat_N1b_ML", "HeliBat_N1c_ML", "Helitron_R25_ML", "Helitron-27_EF", "Helitron-3_EF", "Helitron-33_EF", "Helitron-36_EF", "Helitron-37_EF", "Helitron-39B_EF", "Helitron-52_EF", "Helitron-54_EF", "Helitron-63_EF", "Helitron-8_EF", "Helitron-9_EF"],
	"l1": ["L1MAB_ML", "L1MAB2_ML"],
	"tcmariner": ["MARIN1_ML", "Mariner1_ML", "Mariner3_Ml", "nMar.Tc1_311_ML", "nMar382_Ml", "nMar91_ML", "nMariner-1_EF", "nTc1_135_ML", "nTc1_27_ML", "Tc1_94_ML"],
	"piggybac": ["npiggy_156_ML", "npiggy165_ML", "npiggy2_41_ML", "npiggy259_ML", "piggyBac2_ML"],
	"ves": ["Ves1", "Ves12", "Ves13", "Ves17", "Ves20", "Ves21", "Ves25", "Ves27", "Ves3", "Ves3_ML", "Ves31", "Ves35", "Ves38"]
};

### Extract differential hits that match each species pair and each TE superfamily
### Differential hit = Species A : 0/1 or 1/1 vs Lucifugus : 0/0; or vice versa
### Write all species-pair hits for all elements w/n each superfamily to appropriate output file

# For each species from the list:
for SPECIES in SPECIES_LIST:
	# For each TE superfamily in the TE dictionary:
	for TE_FAMILY in TE_DICT:
		# Make the appropriate output file given the SPECIES, TE superfamily, and mutation rate
		#OUT_FILE = SPECIES + "_" + TE_FAMILY + "_" + MUT_RATE + ".bed"
		#OUTPUT = os.path.join(OUTDIR, OUT_FILE)
		OUT_FILE1 = SPECIES + "_" + TE_FAMILY + "_SPLIT_" + MUT_RATE + ".bed"
		OUTPUT1 = os.path.join(OUTDIR, OUT_FILE1)
		OUT_FILE2 = SPECIES + "_" + TE_FAMILY + "_DEL_" + MUT_RATE + ".bed"
		OUTPUT2 = os.path.join(OUTDIR, OUT_FILE2)
		# Initialize the output list for this species-pair and TE superfamily
		#FAMILY_SP_PAIR_HITS = []
		#Initialize the SPLIT and DELETION output lists for this species-pair and TE superfamily
		FAMILY_SP_PAIR_SPLIT_HITS = []
		FAMILY_SP_PAIR_DEL_HITS = []
		# For each TE subfamily within the TE superfamily:
		for SUBFAMILY in TE_DICT[TE_FAMILY]:
			# Initialize the list for this species-pair and TE subfamily
			#TE_SPECIES_PAIR_HITS = []
			TE_SPECIES_PAIR_SPLIT_HITS = []
			TE_SPECIES_PAIR_DEL_HITS = []
			# Make dataframes for differential hits for the species pair
			# These will correspond to subsets of SPLIT and DELETION hits, respectively
			df1 = df[(df['SVTYPE'] == SUBFAMILY) & (df[SPECIES] != "0/0") & (df['Lucifugus'] == "0/0")]
			df2 = df[(df['SVTYPE'] == SUBFAMILY) & (df[SPECIES] == "0/0") & (df['Lucifugus'] != "0/0")]
			# Subset pf coordinates for non-reference TE insertions (SPLIT hits)
			df_SPLIT = df1[['#CHROM', 'POS', 'END']]
			# Subset of coordinates for reference TE insertions (DELETION hits)
			df_DEL = df2[['#CHROM', 'POS', 'END']]
			# Convert the dataframe of coordinates into a list of lists
			TE_SPECIES_PAIR_SPLIT_HITS = df_SPLIT.values.tolist()
			# Convert the dataframe of coordinates into a list of lists
			TE_SPECIES_PAIR_DEL_HITS = df_DEL.values.tolist()
			
			# Append to the TE superfamily hits coordinate list
			for SUBFAM_HIT in TE_SPECIES_PAIR_SPLIT_HITS:
				FAMILY_SP_PAIR_SPLIT_HITS.append(SUBFAM_HIT)
			for SUBFAM_HIT in TE_SPECIES_PAIR_DEL_HITS:
				FAMILY_SP_PAIR_DEL_HITS.append(SUBFAM_HIT)
		# Once all TE subfamilies for a given superfamily are appended, sort by CHROM and POS
		FAMILY_SP_PAIR_SPLIT_HITS.sort(key=lambda x:(x[0], x[1]))
		FAMILY_SP_PAIR_DEL_HITS.sort(key=lambda x:(x[0], x[1]))
		# Write out each hit within the TE superfamily to the output file
		with open(OUTPUT1, "a+") as f:
			for HIT in FAMILY_SP_PAIR_SPLIT_HITS:
				f.write("{}\t{}\t{}\n".format(HIT[0], HIT[1], HIT[2]))
		with open(OUTPUT2, "a+") as g:
			for HIT in FAMILY_SP_PAIR_DEL_HITS:
				g.write("{}\t{}\t{}\n".format(HIT[0], HIT[1], HIT[2]))
			
			### To simply have a single dataset of all TE insertions, use the following commands
			# Concatenate the two dataframes
			#df3 = pd.concat([df1, df2])
			# Subset the differential hits dataframe to have just the coordinates
			#df4 = df3[['#CHROM', 'POS', 'END']]
			# Convert the dataframe of coordinates into a list of lists
			#TE_SPECIES_PAIR_HITS = df4.values.tolist()
			
			# For each of the hits in the TE subfamily hit list generated:
			# Append to the TE superfamily hits coordinate list
			#for SUBFAM_HIT in TE_SPECIES_PAIR_HITS:
				#FAMILY_SP_PAIR_HITS.append(SUBFAM_HIT)
		# Once all TE subfamilies for a given superfamily are appended, sort by CHROM and POS
		#FAMILY_SP_PAIR_HITS.sort(key=lambda x:(x[0], x[1]))
		# Write out each hit within the TE superfamily to the output file
		#with open(OUTPUT, "a+") as f:
			#for HIT in FAMILY_SP_PAIR_HITS:
				#f.write("{}\t{}\t{}\n".format(HIT[0], HIT[1], HIT[2]))

