import os
import sys
import pandas as pd

# Gathers results of all simulations of a single simulation batch
# (10 replicates for each of n coverage level)
# and places into a single summary file

DIRECTORY = sys.argv[1]
PREFIX = sys.argv[2]

WDIR = "/lustre/scratch/npaulat/ESAT"

COVERAGE = DIRECTORY.split("_")[1]

WDIR = os.path.join(WDIR, DIRECTORY)

os.chdir(WDIR)

OUTPUT = "{}_all_GT_stats.tsv".format(DIRECTORY)

GT_STAT_LIST = {}

for COV in range(1, int(COVERAGE)+1):
	for REP in range(1,11):
		SUBDIR = "{}_{}_{}".format(PREFIX, str(COV), str(REP))
		FILE = "GT_stats.{}_{}_{}.txt".format(PREFIX, str(COV), str(REP))
		STATS = os.path.join(SUBDIR, FILE)
		STATS = os.path.join(WDIR, STATS)
		with open(STATS) as INPUT:
			CONTENT = INPUT.readlines()
		CONTENT = [x.strip() for x in CONTENT if x.strip() != '']
		ITEM = "{}_{}".format(str(COV), str(REP))
		TP_NO_GT = int(CONTENT[4]) - (int(CONTENT[13]) + int(CONTENT[20]))
		FP_NO_GT = int(CONTENT[8]) - (int(CONTENT[15]) + int(CONTENT[22]))
		GT_STAT_LIST[ITEM] = [COV, REP, int(CONTENT[2]), int(CONTENT[4]), int(CONTENT[8]), int(CONTENT[13]), int(CONTENT[20]), TP_NO_GT, int(CONTENT[15]), int(CONTENT[22]), FP_NO_GT]

HEADERS = ['COVERAGE', 'REPLICATE', 'ALL', 'TP', 'FP', 'TP HOM', 'TP HET', 'TP NO GT', 'FP HOM', 'FP HET', 'FP NO GT']

DF = pd.DataFrame.from_dict(GT_STAT_LIST, orient='index', columns=HEADERS)
DF.to_csv(OUTPUT, sep='\t', header=True, index=False)
