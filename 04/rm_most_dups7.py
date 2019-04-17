import sys
import os
import argparse
#import subprocess
#import fnmatch
import itertools
import random
#from Bio import SeqIO
#import copy
import pandas as pd

#Version 7, created 2 Oct 2018 by Nicole S Paulat

##Only for use on MELT Deletion calls, or catenated Deletion+Split calls w/o ASSESS column, which have been modified to include the assumed mLucifugus genotype, and module
#Remove exact duplicate calls, near identical calls, and ambiguous overlapping calls, excluding certain complex cases of near-overlaps or overlaps of TEs of different types or same types with non-close start positions

def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="General removal of most MELT duplicate calls and overlapping calls", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	required = parser.add_argument_group('required arguments')
	#Give input coordinate file of MELT hits (CHROM# + POS)
	parser.add_argument('-i', '--input', help='formatted MELT coordinate list to use as input file in the given directory', required=True)
	#Argument of directory containing the list of formatted MELT hits (need full path)
	parser.add_argument('-d', '--directory', type=str, help='Location of directory of the query scaffold files', required=True)
	#Argument of the output directory (need full path)
	parser.add_argument('-od', '--outdir', type=str, help='Location of directory for the output file', required=True)
	
	args = parser.parse_args()
	COORDS = args.input
	DIR = args.directory
	OUTDIR = args.outdir
	
	#argument sanity checks
	if args.input:
		print('The input coordinate list is ' + COORDS +'.')
	else:
		sys.exit('You must provide the list of MELT coordinates.')
	
	return COORDS, DIR, OUTDIR

COORDS, DIR, OUTDIR = get_args()

#set up access to input file w/ path
INPUT = os.path.join(DIR, COORDS)

#set up access to ouput files w/ path
OUTPUT = os.path.join(OUTDIR, "dupsless_melthits.bed")
OUTPUT2 = os.path.join(OUTDIR, "ambig_melthits.bed")

#open an ouput coordinate file
OUT = open(OUTPUT, "w+")

#open output file for ambiguous MELT calls
CHECK = open(OUTPUT2, "w+")

#cd DIR

#set up function to read in input file into a dictionary
def read_melt(FILENAME):
	my_dict = {}
	for line in open(FILENAME):
		line = line.split()
		CHROM = line[0]
		START = int(line[1])
		STOP = int(line[2])
		TYPE = line[3]
		LENGTH = int(line[4])
		ORIENTATION = line[5]
		GTYPE1 = line[6]
		GTYPE2 = line[7]
		GTYPE3 = line[8]
		GTYPE4 = line[9]
		GTYPE5 = line[10]
		GTYPE6 = line[11]
		GTYPE7 = line[12]
		GTYPE8 = line[13]
		GTYPE9 = line[14]
		GTYPE10 = line[15]
		GTYPE11 = line[16]
		GTYPE12 = line[17]
		MODULE = line[18]
		#if the chromosome is already in the dictionary, add base pos only
		#if the chromosome is not in the dictionary, add entire entry
		try:
			my_dict[line[0]].append((START,STOP,TYPE,LENGTH,ORIENTATION,GTYPE1,GTYPE2,GTYPE3,GTYPE4,GTYPE5,GTYPE6,GTYPE7,GTYPE8,GTYPE9,GTYPE10,GTYPE11,GTYPE12,MODULE))
		except KeyError:
			my_dict[CHROM] = [(START,STOP,TYPE,LENGTH,ORIENTATION,GTYPE1,GTYPE2,GTYPE3,GTYPE4,GTYPE5,GTYPE6,GTYPE7,GTYPE8,GTYPE9,GTYPE10,GTYPE11,GTYPE12,MODULE)]
	return my_dict

#make dictionary from input file
my_dict = read_melt(INPUT)
#my_dict = read_melt("dups_test2.txt")
print('Made first dictionary.')

#initialize output list
output_list = []

#initialize secondary output list - list of calls to manually check with aligned sequence data
confounded_list = []

#make function defining what is within each hit to put on the output list
def melt_call(my_dict,CHROM,INSERTION1):
	return CHROM,my_dict[CHROM][INSERTION1][0],my_dict[CHROM][INSERTION1][1],my_dict[CHROM][INSERTION1][2],my_dict[CHROM][INSERTION1][3],my_dict[CHROM][INSERTION1][4],my_dict[CHROM][INSERTION1][5],my_dict[CHROM][INSERTION1][6],my_dict[CHROM][INSERTION1][7],my_dict[CHROM][INSERTION1][8],my_dict[CHROM][INSERTION1][9],my_dict[CHROM][INSERTION1][10],my_dict[CHROM][INSERTION1][11],my_dict[CHROM][INSERTION1][12],my_dict[CHROM][INSERTION1][13],my_dict[CHROM][INSERTION1][14],my_dict[CHROM][INSERTION1][15],my_dict[CHROM][INSERTION1][16],my_dict[CHROM][INSERTION1][17]

#compare lines to each other
for CHROM in my_dict:
	print('Looking at scaffold ' + CHROM + '.')
	#create lists per each cluster
	#where if same scaffold and start w/n 10 bases or overlap in calls
	#append to "cluster" listdups
	#initialize list of lists of candidate clusters
	INSERTION_CANDIDATES = {x:[] for x in range(len(my_dict[CHROM]))}
	#initialize list of insertions that are close to each other
	ALL_CLOSE = []
	#for every unique combination of calls on a given scaffold
	for INSERTION1, INSERTION2 in itertools.combinations(range(len(my_dict[CHROM])),r=2):
		#define closeness range of w/n 10 bases of an insertion's ends
		a = my_dict[CHROM][INSERTION2][0] - 10
		b = my_dict[CHROM][INSERTION2][1] + 10
		c = my_dict[CHROM][INSERTION1][0] - 10
		d = my_dict[CHROM][INSERTION1][1] + 10
		#if overlapping and/or w/n 10 bases of another insertion, add to "close" list, and list of lists of candidates
		if a <= (my_dict[CHROM][INSERTION1][0] or my_dict[CHROM][INSERTION1][1]) <= b or c <= (my_dict[CHROM][INSERTION2][0] or my_dict[CHROM][INSERTION][1]) <= d:
			INSERTION_CANDIDATES[INSERTION1].append(INSERTION2)
			ALL_CLOSE.append(INSERTION2)
	ALL_CLOSE = set(ALL_CLOSE)
	#print(INSERTION_CANDIDATES)
	#for insertion on "close" list, delete it from the list of lists of candidates (removes calls already in a cluster)
	for INSERTION in ALL_CLOSE:
		del INSERTION_CANDIDATES[INSERTION]
	#for each insertion on the list of candidate cluster lists
	for INSERTION in INSERTION_CANDIDATES:
		#print(INSERTION, INSERTION_CANDIDATES[INSERTION])
		#if insertion has no insertions clustering with it (none w/n "close" range), and not on "close" list, send to output
		if len(INSERTION_CANDIDATES[INSERTION]) == 0:
			if INSERTION not in ALL_CLOSE:
				MELT_CALL = melt_call(my_dict,CHROM,INSERTION)
				output_list.append(MELT_CALL)
		#otherwise, it has a cluster of "close" calls, make further filtering decisions
		else: 
			#recall the actual line data for all insertion calls in a cluster
			INSERTION_CLUSTER = [my_dict[CHROM][INSERTION]] + [my_dict[CHROM][INSERTION2] for INSERTION2 in INSERTION_CANDIDATES[INSERTION]]
			#double-check that there is more than one insertion on a cluster list
			if len(INSERTION_CLUSTER) > 1:
				#make dataframe of the insertion cluster data
				df = pd.DataFrame(INSERTION_CLUSTER)
				#if any hits overlap:
				if all(x[0] <= (z[0] or z[1]) <= x[1] or z[0] <= (x[0] or x[1]) <= z[1] for x,z in itertools.combinations(INSERTION_CLUSTER,r=2)):
				# if orientation is same in all:
					if len(df[4].unique()) == 1:
						#if same genotype per species across all hits:
						if all(len(df[column].unique()) == 1 for column in df.columns[5:16]):
							#print('Meets filter reqs.')
							#if START is same in all or all are w/n 15 bases dist:
							if len(df[0].unique())== 1 or all(abs(x[0]-z[0]) <= 15 for x,z in itertools.combinations(INSERTION_CLUSTER,r=2)):
								#if TE type same in all:
								if len(df[2].unique())== 1:
									#print(df)
									df_subset_MAXLEN = df[df[3] == df[3].max()]
									#print(df_subset_MAXLEN)
									#if multiple hits have max SVLENGTH, pick one at random
									if len(df_subset_MAXLEN) > 1:
										output_list.append([CHROM]+ list(df_subset_MAXLEN.sample(1,axis=0).iloc[0]))
									else:
										output_list.append([CHROM]+ list(df_subset_MAXLEN.iloc[0]))
								#or all are Ves subfamilies:
								elif all(["Ves" in x for x in df[2]]):
									df_subset_MAXLEN = df[df[3] == df[3].max()]
									#if multiple hits have max SVLENGTH, pick one at random
									if len(df_subset_MAXLEN) > 1:
										output_list.append([CHROM]+ list(df_subset_MAXLEN.sample(1,axis=0).iloc[0]))
									else:
										output_list.append([CHROM]+ list(df_subset_MAXLEN.iloc[0]))
								#else if any other TE type combo, add to secondary output for manual review
								else:
									for x in INSERTION_CLUSTER:
										confounded_list.append([CHROM]+ list(x))
							#else if STARTS are more than 15 bases distant, add to secondary output for manual review
							else:
								for x in INSERTION_CLUSTER:
									confounded_list.append([CHROM]+ list(x))
				
				#or if calls are within 10 bases of each other, non-overlapping:
				else:
					continue

# write the final output list to the output file as tab-delimited lines
for line in output_list:
	OUT.write("\t".join((str(x) for x in line)) + '\n')
for line in confounded_list:
	CHECK.write("\t".join((str(x) for x in line)) + '\n')

OUT.close()
CHECK.close()
