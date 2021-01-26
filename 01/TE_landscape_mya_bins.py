import logging
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

# TE Age Landscapes
# Heavily modified from 11 August 2019 script by Jenna R. Grimshaw (jenna.grimshaw@ttu.edu)
# This script will take bedfiles and reformat into separate columns as follows:
# Scaffold, TE, Size, Divergence
# Note: Divergence must be in decimal format, i.e. 1% divergence = 0.01

### Arguments
ap = argparse.ArgumentParser()
ap.add_argument("-b", "--BEDFILE", required=True, help="Please enter your bedfile")
ap.add_argument("-g", "--GENOMESIZE", required=True, type = int, help="Enter size of genome (bp)")
ap.add_argument("-t", "--TAXON", required=True, help="Enter your taxon abbreviation")
ap.add_argument("-f", "--FAMILY", required = False, help = "Choose the TE family you want to plot, e.g. 'Ves', 'MIR', 'L1', 'HAL', 'hAT', 'Helitron', 'piggyBac', 'TcMariner'")
ap.add_argument("-sf", "--SUBFAM_LIST", required = True, help = "List file of the subfamilies you want to plot")
#ap.add_argument("-sf", "--SUPERFAMILY", action = "append", default = [], required = True, help = "Choose the TE superfamily you want to plot, e.g. 'SINE', 'LINE', 'LTR', 'DNA', 'RC'")
#ap.add_argument("-s", "--SUBFAMILY", action = "append", default = [], required = True, help = "Remove false positives that do not contain THIS")
args=vars(ap.parse_args())

##############################################
################ Functions ###################
##############################################

def SUM_OF_SIZE(data, element, min, max):
	size = []
	##Proxy time conversion
	## d = r*t -> d/r = t -> div/(neutral mutation rate) = time (years)
	## div(decimal)/(neutral mut rate per year) = time(years)/1E06 = time(MY)
	## Mammalian neutral mutation rate per year = 2.22 x 10^-09 sub/year
	## So div = 0.01 is ~4.5 MY, div = 0.111 is ~50 MY; div = 0.11 is ~49.55 MY
	rate = 0.00000000222
	mil_year = 1000000
	#CONVERSION = (DIV / rate) / mil_year
	div_min = round(float(min * mil_year * rate),5)
	div_max = round(float(max * mil_year * rate), 5)
	for DIV, SIZE, TE in zip(data.Divergence, data.Size, data.TE):
		if element == TE:
			if DIV > div_min and DIV <= div_max:
				size.append(SIZE)
	return sum(size)
# This function will add up the lengths/sizes 
# of TEs within a min and max Divergence
# and then will return the sum of those lengths

################################################
################################################
################################################

# Enter genome size
genomesize = args["GENOMESIZE"]

# Import RM_bed file
fields = ["Scaffold", "Start", "Stop", "TE", "Length", "Orientation", "Superfamily", "Family", "Divergence", "XX"]
file = pd.read_csv(args["BEDFILE"], sep="\t", names=fields)
print("original file is: ", len(file))

#Import list of TE subfamilies of interested
TELIST = []
NAME = os.path.basename(args["SUBFAM_LIST"]).split(".")[0]
with open(args["SUBFAM_LIST"], "r") as f:
	for line in f:
		line = line.strip()
		TELIST.append(line)

# We are only interested in the FAMILY given
# So create a new dataframe only using TEs from that FAMILY
# Calculate length, drop unneeded columns, reorder columns
if args["FAMILY"]:
	FAMILY = args["FAMILY"]
	if FAMILY == "TcMariner" or "Tc-Mariner" or "TcMar":
		FAMILY == 'TcMar'
		FAM_SUB = file[file.Family.str.contains(FAMILY)]
	else:
		FAM_SUB = file[file.Family == FAMILY]
else:
	FAM_SUB = pd.DataFrame(columns = fields)
	for element in TELIST:
		FAM_SUB = FAM_SUB.append(file[file.TE == element])

FAM_SUB["Size"] = FAM_SUB["Stop"] - FAM_SUB["Start"]
FAM_SUB=FAM_SUB.drop(["Start", "Stop", "Length", "Superfamily", "Orientation", "XX"], axis=1)
FAM_SUB = FAM_SUB[["Scaffold", "Size", "TE","Divergence","Family"]]
FAM_SUB["Divergence"] = FAM_SUB["Divergence"]/100
#print("TE family {} dataframe is: {}".format(FAMILY, len(FAM_SUB)))
print(FAM_SUB.head())

# Remove false positives by only keeping elements that contain "SUBFAMILY"
print("Unique subfamily list:", FAM_SUB["TE"].unique())

# Find unique families
#TELIST = subfamily["TE"].unique()
#TELIST = FAM_SUB["TE"].unique()
print("TE list:", TELIST)

# Create dataframe to be used in loopfile
# I want one column with ranges in increments of millions of years (MY)
# Starting at 1 to 50 by increments of 1
loopfile = pd.DataFrame(columns = TELIST)
my_list = []

for i in range(51):
	my_list.append(i)

loopfile["Range"] = my_list

for TE in TELIST:
	min = 0
	data = []
	for year in loopfile.Range:
		sizesum = SUM_OF_SIZE(FAM_SUB, TE, float(year-1), year)
		data.append(sizesum/genomesize)
	#append TE subfamily data to loopfile
	loopfile[TE] = data
	#print(loopfile)
loopfile.to_csv(args["TAXON"]+"_"+NAME+"_ages.csv", index = False)


for TE in TELIST:
	plt.plot("Range", TE, data = loopfile)
plt.ylabel("Proportion of Genome")
plt.xlabel("Approximate Age (MY)")
##add size parameter to normalize y-axis
#plt.legend()
plt.legend(loc=1)
#plt.show()
plt.savefig(args["TAXON"]+"_"+NAME+".png")
