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

OUTPUT1 = "{}_all_TP_flag_stats.tsv".format(DIRECTORY)
OUTPUT2 = "{}_all_FP_flag_stats.tsv".format(DIRECTORY)
OUTPUT3 = "{}_all_GT_stats.tsv".format(DIRECTORY)

TP_STAT_LIST  = {}

for COV in range(1, int(COVERAGE)+1):
	for REP in range(1,11):
		SUBDIR = "{}_{}_{}".format(PREFIX, str(COV), str(REP))
		FILE = "TP_flag_stats.{}_{}_{}.txt".format(PREFIX, str(COV), str(REP))
		STATS = os.path.join(SUBDIR, FILE)
		STATS = os.path.join(WDIR, STATS)
		with open(STATS) as INPUT:
			CONTENT = INPUT.readlines()
		CONTENT = [x.strip() for x in CONTENT if x.strip() != '']
		ITEM = "{}_{}".format(str(COV), str(REP))
		TP_STAT_LIST[ITEM] = [COV, REP, "TP", int(CONTENT[2]), int(CONTENT[4]), int(CONTENT[6]), int(CONTENT[8]), int(CONTENT[10]), int(CONTENT[12]), int(CONTENT[14]), int(CONTENT[16]), int(CONTENT[18]), int(CONTENT[20]), int(CONTENT[22]), int(CONTENT[24]), int(CONTENT[26]), int(CONTENT[28]), int(CONTENT[30]), int(CONTENT[32]), int(CONTENT[34]), int(CONTENT[36]), int(CONTENT[38]), int(CONTENT[40]), int(CONTENT[42]), int(CONTENT[44]), int(CONTENT[46]), int(CONTENT[48]), int(CONTENT[50]), int(CONTENT[52]), int(CONTENT[54]), int(CONTENT[56]), int(CONTENT[58]), int(CONTENT[60]), int(CONTENT[62]), int(CONTENT[64]), int(CONTENT[66])]

HEADERS1 = ['COVERAGE', 'REPLICATE', 'CALL TYPE', 'Total PASS', 'Total not PASS', 'ac0', 'ac0;s25', 'hDP', 'hDP;ac0', 'hDP;ac0;s25', 'hDP;lc', 'hDP;lc;ac0', 'hDP;lc;ac0;s25', 'hDP;lc;s25', 'hDP;s25', 'lc', 'lc;ac0', 'lc;ac0;s25', 'lc;s25', 'rSD', 'rSD;ac0', 'rSD;ac0;s25', 'rSD;hDP', 'rSD;hDP;ac0', 'rSD;hDP;ac0;s25', 'rSD;hDP;lc', 'rSD;hDP;lc;ac0', 'rSD;hDP;lc;ac0;s25', 'rSD;hDP;lc;s25', 'rSD;hDP;s25', 'rSD;lc', 'rSD;lc;ac0', 'rSD;lc;ac0;s25', 'rSD;lc;s25', 'rSD;s25', 's25']

DF1 = pd.DataFrame.from_dict(TP_STAT_LIST, orient='index', columns=HEADERS1)
DF1.to_csv(OUTPUT1, sep='\t', header=True, index=False)

FP_STAT_LIST  = {}

for COV in range(1, int(COVERAGE)+1):
	for REP in range(1,11):
		SUBDIR = "{}_{}_{}".format(PREFIX, str(COV), str(REP))
		FILE = "FP_flag_stats.{}_{}_{}.txt".format(PREFIX, str(COV), str(REP))
		STATS = os.path.join(SUBDIR, FILE)
		STATS = os.path.join(WDIR, STATS)
		with open(STATS) as INPUT:
			CONTENT = INPUT.readlines()
		CONTENT = [x.strip() for x in CONTENT if x.strip() != '']
		ITEM = "{}_{}".format(str(COV), str(REP))
		FP_STAT_LIST[ITEM] = [COV, REP, "FP", int(CONTENT[2]), int(CONTENT[4]), int(CONTENT[6]), int(CONTENT[8]), int(CONTENT[10]), int(CONTENT[12]), int(CONTENT[14]), int(CONTENT[16]), int(CONTENT[18]), int(CONTENT[20]), int(CONTENT[22]), int(CONTENT[24]), int(CONTENT[26]), int(CONTENT[28]), int(CONTENT[30]), int(CONTENT[32]), int(CONTENT[34]), int(CONTENT[36]), int(CONTENT[38]), int(CONTENT[40]), int(CONTENT[42]), int(CONTENT[44]), int(CONTENT[46]), int(CONTENT[48]), int(CONTENT[50]), int(CONTENT[52]), int(CONTENT[54]), int(CONTENT[56]), int(CONTENT[58]), int(CONTENT[60]), int(CONTENT[62]), int(CONTENT[64]), int(CONTENT[66])]

HEADERS2 = ['COVERAGE', 'REPLICATE', 'CALL TYPE', 'Total PASS', 'Total not PASS', 'ac0', 'ac0;s25', 'hDP', 'hDP;ac0', 'hDP;ac0;s25', 'hDP;lc', 'hDP;lc;ac0', 'hDP;lc;ac0;s25', 'hDP;lc;s25', 'hDP;s25', 'lc', 'lc;ac0', 'lc;ac0;s25', 'lc;s25', 'rSD', 'rSD;ac0', 'rSD;ac0;s25', 'rSD;hDP', 'rSD;hDP;ac0', 'rSD;hDP;ac0;s25', 'rSD;hDP;lc', 'rSD;hDP;lc;ac0', 'rSD;hDP;lc;ac0;s25', 'rSD;hDP;lc;s25', 'rSD;hDP;s25', 'rSD;lc', 'rSD;lc;ac0', 'rSD;lc;ac0;s25', 'rSD;lc;s25', 'rSD;s25', 's25']

DF2 = pd.DataFrame.from_dict(FP_STAT_LIST, orient='index', columns=HEADERS2)
DF2.to_csv(OUTPUT2, sep='\t', header=True, index=False)

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

HEADERS3 = ['COVERAGE', 'REPLICATE', 'ALL', 'TP', 'FP', 'TP HOM', 'TP HET', 'TP NO GT', 'FP HOM', 'FP HET', 'FP NO GT']

DF3 = pd.DataFrame.from_dict(GT_STAT_LIST, orient='index', columns=HEADERS3)
DF3.to_csv(OUTPUT3, sep='\t', header=True, index=False)
