import sys
import os
import argparse
import subprocess
import fnmatch
from Bio import SeqIO

#Compare MELT-DELETION output to its input RepeatMasker file to append orientation
#into the MELT-DELETION hits, and check for inconsistencies in MELT-DEL output

def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="Adding reference genome reference TE orientation to MELT-Deletion hits", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	required = parser.add_argument_group('required arguments')
	#Give input coordinate file of MELT hits (CHROM# + POS + END + SVTYPE)
	parser.add_argument('-i', '--input', help='MELT-Deletion results (in BED format to use as input file in the given RMFILEectory', required=True)
	#Argument of RMFILEectory containing MSA scaffolds (need full path)
	parser.add_argument('-r', '--rmfile', type=str, help='RepeatMasker out.bed file to use as reference input file in the given RMFILEectory', required=True)
	#Argument of RMFILEectory containing MSA scaffolds (need full path)
	parser.add_argument('-od', '--outdir', type=str, help='Location of directory for the output file', required=True)
	
	args = parser.parse_args()
	MDEL = args.input
	RMFILE = args.rmfile
	OUTDIR = args.outdir
	
	return MDEL, RMFILE, OUTDIR
	
MDEL, RMFILE, OUTDIR = get_args()

#argument sanity checks
#if not args.rmfile:
#	sys.exit('You must provide the RepeatMasker file containing the reference TE insertions.')
#if not args.input:
#	sys.exit('You must provide the list of MELT coordinates.')

print('The input MELT Deletion list is ' + MDEL +'.')
print('The input RepeatMasker file is ' + RMFILE + '.')

#cd OUTDIR

#open an ouput file for MELT_DELETION hits with matching RM hits w/n 50 bases
f = open("add_melt_orientation.bed", "w+")

#open output file for MELT-DELETION hits with no corresponding RepeatMasker hits
#either at all, or none within 50 bases ***50 bases is the max cap for matches***
h = open("no_match_meltd.bed", "w+")

#cd RMFILE

def read_melt(FILENAME):
	my_dict = {}
	for line in open(FILENAME):
	#for line in open("test_MDEL.txt"):
		line = line.split()
		CHROM1 = line[0]
		START1 = int(line[1])
		STOP1 = int(line[2])
		TYPE1 = line[3]
		#if the chromosome is already in the dictionary, add base pos only
		#if the chromosome is not in the dictionary, add entire entry
		if line[0] in my_dict:
			my_dict[line[0]].append((START1,STOP1,TYPE1))
			#my_dict[line[0]].append((START1,STOP1,ORIENTATION1))
		else:
			my_dict[CHROM1] = [(START1,STOP1,TYPE1)]
			#my_dict[CHROM1] = [(START1,STOP1,ORIENTATION1)]
	return my_dict		

my_dict = read_melt(MDEL)
#The try/except version of the above:
#	try:
#		my_dict[line[0]].append((START1,STOP1,TYPE1))
#	except KeyError:
#		my_dict[CHROM1] = [(START1,STOP1,TYPE1)]

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
		if line[0] in second_dict:
			second_dict[line[0]].append((START2,STOP2,TYPE2,ORIENTATION2))
			#second_dict[line[0]].append((START1,STOP1,ORIENTATION1))
		else:
			second_dict[CHROM2] = [(START2,STOP2,TYPE2,ORIENTATION2)]
			#second_dict[CHROM2] = [(START2,STOP2,ORIENTATION2)]
	return second_dict
	
second_dict = read_rmout(RMFILE)
#The try/except version of the above:
#	try:
#		second_dict[line[0]].append((START2,STOP2,TYPE2,ORIENTATION2))
#	except KeyError:
#		second_dict[CHROM2] = [(START2,STOP2,TYPE2,ORIENTATION2)]


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
				f.write("{}\t{}\t{}\t{}\t{}\n".format(CHROM,my_dict[CHROM][INSERTION1][0],my_dict[CHROM][INSERTION1][1],my_dict[CHROM][INSERTION1][2],second_dict[CHROM][INSERTION2][3]))
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
				f.write("{}\t{}\t{}\t{}\t{}\n".format(CHROM,my_dict[CHROM][INSERTION1][0],my_dict[CHROM][INSERTION1][1],my_dict[CHROM][INSERTION1][2],BEST_CANDIDATE[3]))
				
			else:
				#INSERTION_CANDIDATES is empty, so write to different file
				h.write("{}\t{}\t{}\t{}\t{}\n".format(CHROM,my_dict[CHROM][INSERTION1][0],my_dict[CHROM][INSERTION1][1],my_dict[CHROM][INSERTION1][2],"NONE"))

#else:
#	print("Unable to read input list. Try again.")
	
#close output file
f.close()
h.close()

#Combine the two files
os.system("cat add_melt_orientation.bed no_match_meltd.bed > all_hits_orient.bed")

def read_catout(FILENAME):
	final_dict = {}
	for line in open(FILENAME):
	#for line in open("all_hits_orient.bed"):
		line = line.split()
		CHROM = line[0]
		START = int(line[1])
		STOP = int(line[2])
		TYPE = line[3]
		ORIENTATION = line[4]
		if line[0] in final_dict:
			final_dict[line[0]].append((START,STOP,TYPE,ORIENTATION))
			#final_dict[line[0]].append((START1,STOP1,ORIENTATION1))
		else:
			final_dict[CHROM] = [(START,STOP,TYPE,ORIENTATION)]
			#final_dict[CHROM2] = [(START2,STOP2,ORIENTATION2)]
	return final_dict

#Read file3 back in as a dictionary
final_dict = read_catout("all_hits_orient.bed")
#open final file for writing
final_file = open("all_hits_orient.bed","w+")

#Loop over the sorted chromosomes and write sorted coordinates
for CHROM in sorted(final_dict.keys()):
	SORTED_COORDS = sorted(final_dict[CHROM], key = lambda x: x[1])
	for INSERTION in SORTED_COORDS:
		#line = line.split()
		#START = int(line[1])
		#STOP = int(line[2])
		#TYPE = line[3]
		#ORIENTATION2 = line[4]
		final_file.write("{}\t{}\t{}\t{}\t{}\n".format(CHROM,INSERTION[0],INSERTION[1],INSERTION[2],INSERTION[3]))

final_file.close()







