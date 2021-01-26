import sys
import os
import argparse
import subprocess
import fnmatch
from Bio import SeqIO

### Version 5, created 27 August 2020 by Nicole S Paulat ###

### Compare MELT-DELETION output to its input RepeatMasker file to append orientation
### into the MELT-DELETION hits, and check for inconsistencies in MELT-DEL output
### Input format is CHROM, POS, END, SVTYPE

def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="Adding reference genome reference TE orientation to MELT-Deletion hits", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	required = parser.add_argument_group('required arguments')
	#Give input coordinate file of MELT hits (CHROM# + POS + END + SVTYPE)
	parser.add_argument('-i', '--input', help='MELT-Deletion results (in BED format to use as input file in the given directory', required=True)
	#Argument of directory containing MSA scaffolds (need full path)
	parser.add_argument('-r', '--rmfile', type=str, help='RepeatMasker out.bed file to use as reference input file in the given directory', required=True)
	#Argument of directory containing the list of formatted MELT hits (need full path)
	parser.add_argument('-d', '--dir', type=str, help='Location of directory of input MELT-Deletion file', default=".")
	#Argument of directory containing MSA scaffolds (need full path)
	parser.add_argument('-od', '--outdir', type=str, help='Location of directory for the output file', default=".")
	
	args = parser.parse_args()
	MDEL = args.input
	RMFILE = args.rmfile
	DIR = args.dir
	OUTDIR = args.outdir
	
	return MDEL, RMFILE, DIR, OUTDIR
	
MDEL, RMFILE, DIR, OUTDIR = get_args()

if DIR == ".":
	DIR = os.getcwd()
if OUTDIR == ".":
	OUTDIR = os.getcwd()

#argument sanity checks
#if not args.rmfile:
#	sys.exit('You must provide the RepeatMasker file containing the reference TE insertions.')
#if not args.input:
#	sys.exit('You must provide the list of MELT coordinates.')

print('The input MELT Deletion list is ' + MDEL +'.')
print('The input RepeatMasker file is ' + RMFILE + '.')

#set up access to input file w/ path
INPUT = os.path.join(DIR, MDEL)

#set up access to ouput files w/ path
BASENAME = os.path.basename(MDEL).split(".")[0]
BASE_OUT = os.path.join(OUTDIR, BASENAME)
OUTPUT = BASE_OUT + "_add_melt_orientation.bed"
OUTPUT2 = BASE_OUT + "_no_match_meltd.bed"

#open an ouput file for MELT_DELETION hits with matching RM hits w/n 50 bases
f = open(OUTPUT, "w+")

#open output file for MELT-DELETION hits with no corresponding RepeatMasker hits
#either at all, or none within 50 bases ***50 bases is the max cap for matches***
h = open(OUTPUT2, "w+")

# def read_melt(FILENAME):
	# my_dict = {}
	# for line in open(FILENAME):
	# #for line in open("test_MDEL.txt"):
		# line = line.split()
		# CHROM1 = line[0]
		# START1 = int(line[1])
		# STOP1 = int(line[2])
		# TYPE1 = line[3]
		# #if the chromosome is already in the dictionary, add base pos only
		# #if the chromosome is not in the dictionary, add entire entry
		# #if line[0] in my_dict:
		# try:
			# my_dict[line[0]].append((START1,STOP1,TYPE1))
		# #else:
		# except KeyError:
			# my_dict[CHROM1] = [(START1,STOP1,TYPE1)]
	# return my_dict

def read_melt(FILENAME):
	my_dict = {}
	for line in open(FILENAME):
		line = line.split()
		CHROM = line[0]
		START = int(line[1])
		STOP = int(line[2])
		TYPE = line[3]
		LENGTH = int(line[4])
		GTYPE1 = line[5]
		GTYPE2 = line[6]
		GTYPE3 = line[7]
		GTYPE4 = line[8]
		GTYPE5 = line[9]
		GTYPE6 = line[10]
		GTYPE7 = line[11]
		GTYPE8 = line[12]
		GTYPE9 = line[13]
		GTYPE10 = line[14]
		GTYPE11 = line[15]
		MODULE = line[16]
		#if the chromosome is already in the dictionary, add base pos only
		#if the chromosome is not in the dictionary, add entire entry
		try:
			my_dict[line[0]].append((START,STOP,TYPE,LENGTH,GTYPE1,GTYPE2,GTYPE3,GTYPE4,GTYPE5,GTYPE6,GTYPE7,GTYPE8,GTYPE9,GTYPE10,GTYPE11,MODULE))
		except KeyError:
			my_dict[CHROM] = [(START,STOP,TYPE,LENGTH,GTYPE1,GTYPE2,GTYPE3,GTYPE4,GTYPE5,GTYPE6,GTYPE7,GTYPE8,GTYPE9,GTYPE10,GTYPE11,MODULE)]
	return my_dict

my_dict = read_melt(INPUT)

print('Made first dictionary.')
		
def read_rmout(FILENAME):
	second_dict = {}
	for line in open(FILENAME):
	#for line in open("test_rmfile.txt"):
		line = line.split()
		CHROM2 = line[0]
		START2 = int(line[1])
		STOP2 = int(line[2])
		TYPE2 = line[3]
		ORIENTATION2 = line[5]
		#if line[0] in second_dict:
		try:
			second_dict[line[0]].append((START2,STOP2,TYPE2,ORIENTATION2))
		#else:
		except KeyError:
			second_dict[CHROM2] = [(START2,STOP2,TYPE2,ORIENTATION2)]
	return second_dict
	
second_dict = read_rmout(RMFILE)

print('Made second dictionary.')

#for every key in the dictionary, compare to every line in RMFILE
#need sorted MDEL list
#edit to list
for CHROM in my_dict:
	print('Looking at scaffold ' + CHROM + '.')
	if CHROM in second_dict:
		print('Scaffold is in RepeatMasker file.')
		#print(my_dict[CHROM])
		#print(second_dict[CHROM])
		for INSERTION1 in range(len(my_dict[CHROM])):
			INSERTION_CANDIDATES = []
			for INSERTION2 in range(len(second_dict[CHROM])):
				if my_dict[CHROM][INSERTION1][2] == second_dict[CHROM][INSERTION2][2]:
					if abs(my_dict[CHROM][INSERTION1][0] - second_dict[CHROM][INSERTION2][0]) <= 50:
						#If insertion meets the criteria, add the insertion (tuple) to a list
						INSERTION_CANDIDATES.append(second_dict[CHROM][INSERTION2])

			if len(INSERTION_CANDIDATES) == 1:
				#Only one candidate, write it out
				#f.write("{}\t{}\t{}\t{}\t{}\n".format(CHROM,my_dict[CHROM][INSERTION1][0],my_dict[CHROM][INSERTION1][1],my_dict[CHROM][INSERTION1][2],second_dict[CHROM][INSERTION2][3]))
				f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(CHROM,my_dict[CHROM][INSERTION1][0],my_dict[CHROM][INSERTION1][1],my_dict[CHROM][INSERTION1][2], my_dict[CHROM][INSERTION1][3], second_dict[CHROM][INSERTION2][3], my_dict[CHROM][INSERTION1][4], my_dict[CHROM][INSERTION1][5], my_dict[CHROM][INSERTION1][6], my_dict[CHROM][INSERTION1][7], my_dict[CHROM][INSERTION1][8], my_dict[CHROM][INSERTION1][9], my_dict[CHROM][INSERTION1][10], my_dict[CHROM][INSERTION1][11], my_dict[CHROM][INSERTION1][12], my_dict[CHROM][INSERTION1][13], my_dict[CHROM][INSERTION1][14], my_dict[CHROM][INSERTION1][15]))
						#sys.stdout.write("{}\t{}\t{}\t{}\t{}\n".format(CHROM,my_dict[CHROM][INSERTION1][0],my_dict[CHROM][INSERTION1][1],my_dict[CHROM][INSERTION1][2],second_dict[CHROM][INSERTION2][3]))
			elif len(INSERTION_CANDIDATES) > 1:
				#Code for choosing among insertion candidates
				INSERTION_START = my_dict[CHROM][INSERTION1][0]
				# list comprehension to calculate difference between MELT start and RM start for all candidates
				START_DIFFS = [abs(INSERTION_START - I[0]) for I in INSERTION_CANDIDATES]
				# This is a way to get the index of the minimum in a list
				# NOTE: It will always return the first instance of the minimum
				MIN_START_DIFF = START_DIFFS.index(min(START_DIFFS))
				BEST_CANDIDATE = INSERTION_CANDIDATES[MIN_START_DIFF]
				#f.write("{}\t{}\t{}\t{}\t{}\n".format(CHROM,my_dict[CHROM][INSERTION1][0],my_dict[CHROM][INSERTION1][1],my_dict[CHROM][INSERTION1][2],BEST_CANDIDATE[3]))
				f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(CHROM,my_dict[CHROM][INSERTION1][0],my_dict[CHROM][INSERTION1][1],my_dict[CHROM][INSERTION1][2], my_dict[CHROM][INSERTION1][3], BEST_CANDIDATE[3], my_dict[CHROM][INSERTION1][4], my_dict[CHROM][INSERTION1][5], my_dict[CHROM][INSERTION1][6], my_dict[CHROM][INSERTION1][7], my_dict[CHROM][INSERTION1][8], my_dict[CHROM][INSERTION1][9], my_dict[CHROM][INSERTION1][10], my_dict[CHROM][INSERTION1][11], my_dict[CHROM][INSERTION1][12], my_dict[CHROM][INSERTION1][13], my_dict[CHROM][INSERTION1][14], my_dict[CHROM][INSERTION1][15]))
			else:
				#INSERTION_CANDIDATES is empty, so write to different file
				#h.write("{}\t{}\t{}\t{}\t{}\n".format(CHROM,my_dict[CHROM][INSERTION1][0],my_dict[CHROM][INSERTION1][1],my_dict[CHROM][INSERTION1][2],"NONE"))
				h.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(CHROM,my_dict[CHROM][INSERTION1][0],my_dict[CHROM][INSERTION1][1],my_dict[CHROM][INSERTION1][2], my_dict[CHROM][INSERTION1][3], "NONE", my_dict[CHROM][INSERTION1][4], my_dict[CHROM][INSERTION1][5], my_dict[CHROM][INSERTION1][6], my_dict[CHROM][INSERTION1][7], my_dict[CHROM][INSERTION1][8], my_dict[CHROM][INSERTION1][9], my_dict[CHROM][INSERTION1][10], my_dict[CHROM][INSERTION1][11], my_dict[CHROM][INSERTION1][12], my_dict[CHROM][INSERTION1][13], my_dict[CHROM][INSERTION1][14], my_dict[CHROM][INSERTION1][15]))

#else:
#	print("Unable to read input list. Try again.")
	
#close output files
f.close()
h.close()

#Combine the two files
CAT_FILE = BASE_OUT + "_all_hits_orient.bed"
CONCATENATE = "cat " + OUTPUT + " " + OUTPUT2 + " > " + CAT_FILE
os.system(CONCATENATE)

def read_catout(FILENAME):
	final_dict = {}
	#for line in open("all_hits_orient.bed"):
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
		MODULE = line[17]
		#if the chromosome is already in the dictionary, add base pos only
		#if the chromosome is not in the dictionary, add entire entry
		try:
			final_dict[line[0]].append((START,STOP,TYPE,LENGTH,ORIENTATION,GTYPE1,GTYPE2,GTYPE3,GTYPE4,GTYPE5,GTYPE6,GTYPE7,GTYPE8,GTYPE9,GTYPE10,GTYPE11,MODULE))
		except KeyError:
			final_dict[CHROM] = [(START,STOP,TYPE,LENGTH,ORIENTATION,GTYPE1,GTYPE2,GTYPE3,GTYPE4,GTYPE5,GTYPE6,GTYPE7,GTYPE8,GTYPE9,GTYPE10,GTYPE11,MODULE)]
	return final_dict
		# line = line.split()
		# CHROM = line[0]
		# START = int(line[1])
		# STOP = int(line[2])
		# TYPE = line[3]
		# ORIENTATION = line[4]
		#if line[0] in final_dict:
		# try:
			# final_dict[line[0]].append((START,STOP,TYPE,ORIENTATION))
		# #else:
		# except KeyError:
			# final_dict[CHROM] = [(START,STOP,TYPE,ORIENTATION)]
	# return final_dict

#Read file3 back in as a dictionary
final_dict = read_catout(CAT_FILE)
#open final file for writing
final_file = open(CAT_FILE,"w+")

#Loop over the sorted chromosomes and write sorted coordinates
for CHROM in sorted(final_dict.keys()):
	SORTED_COORDS = sorted(final_dict[CHROM], key = lambda x: x[1])
	for INSERTION in SORTED_COORDS:
		#final_file.write("{}\t{}\t{}\t{}\t{}\n".format(CHROM,INSERTION[0],INSERTION[1],INSERTION[2],INSERTION[3]))
		final_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(CHROM,INSERTION[0],INSERTION[1],INSERTION[2],INSERTION[3], INSERTION[4], INSERTION[5], INSERTION[6], INSERTION[7], INSERTION[8], INSERTION[9], INSERTION[10], INSERTION[11], INSERTION[12], INSERTION[13], INSERTION[14], INSERTION[15], INSERTION[16]))

final_file.close()
