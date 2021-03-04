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

STAT_LIST  = {}

for COV in range(1, int(COVERAGE)+1):
	for REP in range(1,11):
		SUBDIR = "{}_{}_{}".format(PREFIX, str(COV), str(REP))
		FILE = "stats.{}_{}_{}.txt".format(PREFIX, str(COV), str(REP))
		STATS = os.path.join(SUBDIR, FILE)
		STATS = os.path.join(WDIR, STATS)
		with open(STATS) as INPUT:
			CONTENT = INPUT.readlines()
		CONTENT = [x.strip() for x in CONTENT if x.strip() != '']
		ITEM = "{}_{}".format(str(COV), str(REP))
		STAT_LIST[ITEM] = [COV, REP, int(CONTENT[7]), int(CONTENT[11]), int(CONTENT[9]), int(CONTENT[16]), int(CONTENT[20]), int(CONTENT[18]), int(CONTENT[25]), int(CONTENT[29]), int(CONTENT[27]), int(CONTENT[34]), int(CONTENT[38]), int(CONTENT[36])]

OUTPUT = "{}_all_stats.tsv".format(DIRECTORY)

HEADERS = ['COVERAGE', 'REPLICATE', 'ALL_TP', 'ALL_FP', 'ALL_FN', 'PASS_TP', 'PASS_FP', 'PASS_FN', 'PASS_HOM_TP', 'PASS_HOM_FP', 'PASS_HOM_FN', 'PASS_HET_TP', 'PASS_HET_FP', 'PASS_HET_FN']

DF = pd.DataFrame.from_dict(STAT_LIST, orient='index', columns=HEADERS)
DF.to_csv(OUTPUT, sep='\t', header=True, index=False)
