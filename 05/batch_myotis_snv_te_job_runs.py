import sys
import os
import argparse
#import itertools
#import subprocess

#Arguments: #species#, vcf type (snp or indel), TE type, mutation rate, outdir, queue

def get_args():
	parser = argparse.ArgumentParser(description="Batch script generator for SNP-TE R analysis runs given MELT-derived coordinate bed files for a specific max TE mutation rate", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	#required = parser.add_argument_group('required arguments')
	parser.add_argument('-r', '--refindex', type=str, help='Name of the reference genome\'s index file (.fai)', required=True)
	parser.add_argument('-m', '--rate', type=int, help="Maximum mutation rate of the MELT run the coordinates are taken from", required=True)
	parser.add_argument('-np', '--proc', type=int, help='Number of cores to use for multithreaded applications', required=True)
	parser.add_argument('-d', '--dir', type=str, help="Path of the directory for the input MELT hit bed files, need full path", required=True)
	parser.add_argument('-wd', '--wdir', type=str, help="Full path of directory containing all of the VCFs, FAI, and Rscript template files", required=True)
	parser.add_argument('-od', '--outdir', type=str, help='Location of directory for output files', default=".")
	parser.add_argument('-q', '--queue', type=str, help='quanah or hrothgar', required=True)

	args = parser.parse_args()
	REF = args.refindex
	RATE = args.rate
	PROC = args.proc
	DIR = args.dir
	WDIR = args.wdir
	OUTDIR = args.outdir
	QUEUE = args.queue
	
	return REF, RATE, PROC, DIR, WDIR, OUTDIR, QUEUE
	
REF, RATE, PROC, DIR, WDIR, OUTDIR, QUEUE = get_args()
#VMEM = 6

print('The reference genome index file is ' + REF + '.')
print('The max mutation rate of these MELT coords is ' + str(RATE) + ' mutations in 100 bp.')
print('Use ' + str(PROC) + ' processors.')
print('The MELT hits BED file directory is ' + DIR + '.')
print('The input VCFs and FAI files are in ' + WDIR + '.')
print('The output directory is ' + OUTDIR + '.')
print('The queue is ' + QUEUE + '.') 

if OUTDIR == '.':
	OUTDIR = os.getcwd()

MUT_RATE = 'mut' + str(RATE)
MUT_RATE_DIR = os.path.join(OUTDIR, MUT_RATE)
if not os.path.exists(MUT_RATE_DIR):
	print(MUT_RATE_DIR + ' does not exist. Make it.')
	os.mkdir(MUT_RATE_DIR)
else:
	print(MUT_RATE_DIR + ' exists.')

INPUT_BED_NAME = os.path.abspath(DIR).split('/')[-1]
OUT_SUBDIR = os.path.join(MUT_RATE_DIR, INPUT_BED_NAME)
if os.path.exists(OUT_SUBDIR):
	print("Input subdirectory {} exists.".format(OUT_SUBDIR))
else:
	os.mkdir(OUT_SUBDIR)
os.chdir(OUT_SUBDIR)

INDEX_FILE = os.path.basename(REF)
REF_INDEX = os.path.join(WDIR, INDEX_FILE)

TEMPLATE_RSCRIPT = "snp_te_analysis_template7.r"
TEMPLATE_PATH = os.path.join(WDIR, TEMPLATE_RSCRIPT)

SCRIPT_DIR = "scripts"
SCRIPT_SUBDIR = os.path.join(OUT_SUBDIR, SCRIPT_DIR)
if os.path.exists(SCRIPT_SUBDIR):
	print("Individual submission scripts subdirectory {} exists.".format(SCRIPT_SUBDIR))
else:
	os.mkdir(SCRIPT_SUBDIR)

SPECIES_LIST = ['Austroriparius', 'Brandtii', 'Ciliolabrum', 'Davidii', 'Occultus', 'Sept_TTU', 'Thysanodes', 'Velifer', 'Vivesi', 'Yumanensis']

#with open("/lustre/scratch/npaulat/MELTv2.1.5/references/te_list.txt", "r") as d:
#	TES = d.read().split(" ")
TE_FAMILY_LIST = ['HAL', 'LINE', 'SINE', 'hAT', 'PiggyBac', 'TcMariner', 'Helitron']

VARIANT_LIST = ['snps', 'indels']

if "cat" in INPUT_BED_NAME:
	IN_GROUP_LIST = ['all', 'SPLIT', 'DEL']
else:
	IN_GROUP_LIST = ['all']

for VTYPE in VARIANT_LIST:
	MUT_RATE_VTYPE_DIR = os.path.join(OUT_SUBDIR, VTYPE)
	if not os.path.exists(MUT_RATE_VTYPE_DIR):
		print(MUT_RATE_VTYPE_DIR + ' does not exist. Make it.')
		os.mkdir(MUT_RATE_VTYPE_DIR)
	else:
		print(MUT_RATE_VTYPE_DIR + ' exists.')
	
	JOBSSUBMISSON = 'jobs_submission_all_' + QUEUE + '_' + INPUT_BED_NAME + '_' + VTYPE +  '.sh'
	JOBSUB_FILE = os.path.join(SCRIPT_SUBDIR, JOBSSUBMISSON)
	
	for TE_TYPE in TE_FAMILY_LIST:
		for SPECIES in SPECIES_LIST:
			for IN_GROUP in IN_GROUP_LIST:
				QSUBFILENAME = 'qsub_' + QUEUE + '_' + INPUT_BED_NAME + '_' + SPECIES + '_' + IN_GROUP + '_' + VTYPE + '_' + TE_TYPE + '.sh'
				QSUB_FILE = os.path.join(SCRIPT_SUBDIR, QSUBFILENAME)
				RSCRIPTNAME = INPUT_BED_NAME + '_' + SPECIES + '_' + IN_GROUP + '_' + VTYPE + '_' + TE_TYPE + '.r'
				RSCRIPT_FILE = os.path.join(OUT_SUBDIR, RSCRIPTNAME)
				BED_FILE = SPECIES + '_' + TE_TYPE + '_' + IN_GROUP + '_' + MUT_RATE + '.bed'
				BED_PATH = os.path.join(DIR, BED_FILE)
				VCF_FILE = SPECIES + '_' + VTYPE + '.vcf'
				VCF_PATH = os.path.join(WDIR, VCF_FILE)
				
				#For Hrothgar queue.
				if QUEUE == 'hrothgar':
					with open(QSUB_FILE, 'w+') as f:
						f.write('#!/bin/sh' + '\n')
						f.write('#$ -V' + '\n')
						f.write('#$ -cwd' + '\n')
						f.write('#$ -S /bin/bash' + '\n')
						f.write('#$ -N ' + VTYPE + '_' + TE_TYPE + '_' + MUT_RATE + '\n')
						f.write('#$ -o $JOB_NAME.o$JOB_ID' + '\n')
						f.write('#$ -e $JOB_NAME.e$JOB_ID' + '\n')
						f.write('#$ -q Chewie' + '\n')
						f.write('#$ -pe fill ' + str(PROC) + '\n')
						f.write('#$ -P communitycluster' + '\n')
						#f.write('#$ -l h_vmem=' + str(VMEM) + 'G' + '\n')
						f.write('\n')
						f.write('#The reference genome index is ' + REF + '.\n')
						f.write('#The max mutation rate from MELT is ' + str(RATE) + ' mutations in 100 bp.\n')
						#f.write('#Use ' + str(PROC) + ' processors and ' + str(VMEM) + 'G memory.\n')
						f.write('#Use ' + str(PROC) + ' processors.\n')
						f.write('#The output directory is ' + MUT_RATE_VTYPE_DIR + '/del.\n')
						f.write('#The queue is ' + QUEUE + '.\n') 
						f.write('\n')
						f.write('module load intel/18.0.3.222 impi/2018.3.222 R-3.6.1/3.6.1-intel\n')
						f.write('\n')
						f.write('#Change all variables in R script\n')
						f.write('VCF_FILE="' + VCF_PATH + '"\n')
						f.write('BED_FILE="' + BED_PATH + '"\n')
						f.write('DATASET="' + INPUT_BED_NAME + '"\n')
						f.write('OUTDIR="' + MUT_RATE_VTYPE_DIR + '"\n')
						f.write('SPECIES="' + SPECIES + '"\n')
						f.write('IN_TYPE="' + IN_GROUP + '"\n')
						f.write('REFINDEX="' + REF_INDEX + '"\n')
						f.write('TE_TYPE="' + TE_TYPE + '"\n')
						f.write('VTYPE="' + VTYPE + '"\n')
						if VTYPE == 'snps':
							f.write('VARIANT="SNPs"' + '\n')
						elif VTYPE == 'indels':
							f.write('VARIANT="Indels"' + '\n')
						f.write('sed "s/<SPECIES>/$SPECIES/g" ' + TEMPLATE_PATH + ' > ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<VCF_FILE>|$VCF_FILE|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<BED_FILE>|$BED_FILE|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<DATASET>|$DATASET|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<OUTDIR>|$OUTDIR|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<IN_TYPE>|$IN_TYPE|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<REFINDEX>|$REFINDEX|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<TE_TYPE>|$TE_TYPE|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<VTYPE>|$VTYPE|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<VARIANT>|$VARIANT|g" ' + RSCRIPT_FILE + '\n')
						f.write('\n')
						f.write('Rscript ' + RSCRIPT_FILE + '\n')
						with open(JOBSUB_FILE, 'a+') as g:
							g.write('qsub ' + QSUB_FILE + '\n')

				#For Quanah queue.
				elif QUEUE == 'quanah':
					with open(QSUB_FILE, 'w+') as f:
						f.write('#!/bin/sh' + '\n')
						f.write('#$ -V' + '\n')
						f.write('#$ -cwd' + '\n')
						f.write('#$ -S /bin/bash' + '\n')
						f.write('#$ -N ' + VTYPE + '_' + TE_TYPE + '_' + MUT_RATE + '\n')
						f.write('#$ -o $JOB_NAME.o$JOB_ID' + '\n')
						f.write('#$ -e $JOB_NAME.e$JOB_ID' + '\n')
						f.write('#$ -q omni' + '\n')
						f.write('#$ -pe sm ' + str(PROC) + '\n')
						f.write('#$ -P quanah' + '\n')
						#f.write('#$ -l h_vmem=' + str(VMEM) + 'G' + '\n')
						f.write('\n')
						f.write('#The reference genome index is ' + REF + '.\n')
						f.write('#The max mutation rate from MELT is ' + str(RATE) + ' mutations in 100 bp.\n')
						#f.write('#Use ' + str(PROC) + ' processors and ' + str(VMEM) + 'G memory.\n')
						f.write('#Use ' + str(PROC) + ' processors.\n')
						f.write('#The output directory is ' + MUT_RATE_VTYPE_DIR + '/del.\n')
						f.write('#The queue is ' + QUEUE + '.\n') 
						f.write('\n')
						f.write('module load intel/18.0.3.222 impi/2018.3.222 R/3.6.1-intel\n')
						f.write('\n')
						f.write('#Change all variables in R script\n')
						f.write('VCF_FILE="' + VCF_PATH + '"\n')
						f.write('BED_FILE="' + BED_PATH + '"\n')
						f.write('DATASET="' + INPUT_BED_NAME + '"\n')
						f.write('OUTDIR="' + MUT_RATE_VTYPE_DIR + '"\n')
						f.write('SPECIES="' + SPECIES + '"\n')
						f.write('IN_TYPE="' + IN_GROUP + '"\n')
						f.write('REFINDEX="' + REF_INDEX + '"\n')
						f.write('TE_TYPE="' + TE_TYPE + '"\n')
						f.write('VTYPE="' + VTYPE + '"\n')
						if VTYPE == 'snps':
							f.write('VARIANT="SNPs"' + '\n')
						elif VTYPE == 'indels':
							f.write('VARIANT="Indels"' + '\n')
						f.write('sed "s/<SPECIES>/$SPECIES/g" ' + TEMPLATE_PATH + ' > ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<VCF_FILE>|$VCF_FILE|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<BED_FILE>|$BED_FILE|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<DATASET>|$DATASET|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<OUTDIR>|$OUTDIR|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<IN_TYPE>|$IN_TYPE|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<REFINDEX>|$REFINDEX|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<TE_TYPE>|$TE_TYPE|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<VTYPE>|$VTYPE|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<VARIANT>|$VARIANT|g" ' + RSCRIPT_FILE + '\n')
						f.write('\n')
						f.write('Rscript ' + RSCRIPT_FILE + '\n')
						with open(JOBSUB_FILE, 'a+') as g:
							g.write('qsub ' + QSUB_FILE + '\n')
				else:
					print('Bad queue choice. Your only choices are hrothgar and quanah.')
