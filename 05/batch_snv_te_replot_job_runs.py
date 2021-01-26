import sys
import os
import argparse

# For replotting SNV-TE graphs that were previously made using batch_myotis_snv_te_job_runs.py

#Arguments: Mutation rate, # processors, working directory, output directory, queue

def get_args():
	parser = argparse.ArgumentParser(description="Batch script generator for SNP-TE R analysis runs given MELT-derived coordinate bed files for a specific max TE mutation rate", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	#required = parser.add_argument_group('required arguments')
	parser.add_argument('-m', '--rate', type=int, help="Maximum mutation rate of the MELT run the coordinates are taken from", required=True)
	parser.add_argument('-np', '--proc', type=int, help='Number of cores to use for multithreaded applications', required=True)
	parser.add_argument('-wd', '--wdir', type=str, help="Full path of directory containing the Rscript template files", required=True)
	parser.add_argument('-od', '--outdir', type=str, help='Location of directory for output files', default=".")
	parser.add_argument('-q', '--queue', type=str, help='quanah or hrothgar', required=True)

	args = parser.parse_args()
	RATE = args.rate
	PROC = args.proc
	WDIR = args.wdir
	OUTDIR = args.outdir
	QUEUE = args.queue
	
	return RATE, PROC, WDIR, OUTDIR, QUEUE
	
RATE, PROC, WDIR, OUTDIR, QUEUE = get_args()

print('The max mutation rate of these MELT coords is ' + str(RATE) + ' mutations in 100 bp.')
print('Use ' + str(PROC) + ' processors.')
print('The input RData files are in ' + WDIR + '.')
print('The output directory is ' + OUTDIR + '.')
print('The queue is ' + QUEUE + '.') 

if OUTDIR == '.':
	OUTDIR = os.getcwd()

MUT_RATE = 'mut' + str(RATE)

INPUT_NAME = "cat_final"
INPUT_SUBDIR = "mut10/mut10_cat_final"
OUT_SUBDIR = os.path.join(OUTDIR, INPUT_SUBDIR)
os.chdir(OUT_SUBDIR)

TEMPLATE_RSCRIPT = "snp_te_replot_template_snps.r"
#TEMPLATE_RSCRIPT = "snp_te_replot_template_indels.r"
TEMPLATE_PATH = os.path.join(WDIR, TEMPLATE_RSCRIPT)

SCRIPT_DIR = "scripts"
SCRIPT_SUBDIR = os.path.join(OUT_SUBDIR, SCRIPT_DIR)
if os.path.exists(SCRIPT_SUBDIR):
	print("Individual submission scripts subdirectory {} exists.".format(SCRIPT_SUBDIR))
#else:
#	os.mkdir(SCRIPT_SUBDIR)

SPECIES_LIST = ['Austroriparius', 'Brandtii', 'Ciliolabrum', 'Davidii', 'Occultus', 'Sept_TTU', 'Thysanodes', 'Velifer', 'Vivesi', 'Yumanensis']

TE_FAMILY_LIST = TE_FAMILY_LIST = ['HAL', 'LINE', 'SINE', 'hAT', 'PiggyBac', 'TcMariner', 'Helitron']

VARIANT_LIST = ['snps', 'indels']

if "cat" in INPUT_NAME:
#	IN_GROUP_LIST = ['all', 'SPLIT', 'DEL']
#else:
	IN_GROUP_LIST = ['all']

for VTYPE in VARIANT_LIST:
	MUT_RATE_VTYPE_DIR = os.path.join(OUT_SUBDIR, VTYPE)
	
	JOBSSUBMISSON = 'jobs_submission_all_' + QUEUE + '_' + INPUT_NAME + '_' + VTYPE + '_replot.sh'
	JOBSUB_FILE = os.path.join(SCRIPT_SUBDIR, JOBSSUBMISSON)
	
	for TE_TYPE in TE_FAMILY_LIST:
		for SPECIES in SPECIES_LIST:
			for IN_GROUP in IN_GROUP_LIST:
				QSUBFILENAME = 'qsub_' + QUEUE + '_' + INPUT_NAME + '_' + SPECIES + '_' + IN_GROUP + '_' + VTYPE + '_' + TE_TYPE + '_replot.sh'
				QSUB_FILE = os.path.join(SCRIPT_SUBDIR, QSUBFILENAME)
				RSCRIPTNAME = INPUT_NAME + '_' + SPECIES + '_' + IN_GROUP + '_' + VTYPE + '_' + TE_TYPE + '_replot.r'
				RSCRIPT_FILE = os.path.join(OUT_SUBDIR, RSCRIPTNAME)
				RDATA_FILE = MUT_RATE + '_' + INPUT_NAME + '_' + SPECIES + '_' + TE_TYPE + '_' + IN_GROUP + '_' + VTYPE + '.RData'
				RDATA_PATH = os.path.join(MUT_RATE_VTYPE_DIR, RDATA_FILE)
				
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
						f.write('#$ -q Yoda' + '\n')
						f.write('#$ -pe sm ' + str(PROC) + '\n')
						f.write('#$ -P communitycluster' + '\n')
						#f.write('#$ -l h_vmem=' + str(VMEM) + 'G' + '\n')
						f.write('\n')
						f.write('#The max mutation rate from MELT is ' + str(RATE) + ' mutations in 100 bp.\n')
						#f.write('#Use ' + str(PROC) + ' processors and ' + str(VMEM) + 'G memory.\n')
						f.write('#Use ' + str(PROC) + ' processors.\n')
						f.write('#The output directory is ' + MUT_RATE_VTYPE_DIR + '/.\n')
						f.write('#The queue is ' + QUEUE + '.\n') 
						f.write('\n')
						f.write('module load intel/18.0.3.222 impi/2018.3.222 R-3.6.1/3.6.1-intel\n')
						f.write('\n')
						f.write('#Change all variables in R script\n')
						f.write('DATASET="' + MUT_RATE + "_" + INPUT_NAME + '"\n')
						f.write('OUTDIR="' + MUT_RATE_VTYPE_DIR + '"\n')
						f.write('SPECIES="' + SPECIES + '"\n')
						f.write('IN_TYPE="' + IN_GROUP + '"\n')
						f.write('TE_TYPE="' + TE_TYPE + '"\n')
						f.write('VTYPE="' + VTYPE + '"\n')
						if VTYPE == 'snps':
							f.write('VARIANT="SNPs"' + '\n')
						elif VTYPE == 'indels':
							f.write('VARIANT="Indels"' + '\n')
						f.write('sed "s/<SPECIES>/$SPECIES/g" ' + TEMPLATE_PATH + ' > ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<DATASET>|$DATASET|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<OUTDIR>|$OUTDIR|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<IN_TYPE>|$IN_TYPE|g" ' + RSCRIPT_FILE + '\n')
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
						f.write('#The max mutation rate from MELT is ' + str(RATE) + ' mutations in 100 bp.\n')
						#f.write('#Use ' + str(PROC) + ' processors and ' + str(VMEM) + 'G memory.\n')
						f.write('#Use ' + str(PROC) + ' processors.\n')
						f.write('#The output directory is ' + MUT_RATE_VTYPE_DIR + '/del.\n')
						f.write('#The queue is ' + QUEUE + '.\n') 
						f.write('\n')
						f.write('module load intel/18.0.3.222 impi/2018.3.222 R/3.5.0-intel\n')
						f.write('\n')
						f.write('#Change all variables in R script\n')
						f.write('DATASET="' + MUT_RATE + "_" + INPUT_NAME + '"\n')
						f.write('OUTDIR="' + MUT_RATE_VTYPE_DIR + '"\n')
						f.write('SPECIES="' + SPECIES + '"\n')
						f.write('IN_TYPE="' + IN_GROUP + '"\n')
						f.write('TE_TYPE="' + TE_TYPE + '"\n')
						f.write('VTYPE="' + VTYPE + '"\n')
						if VTYPE == 'snps':
							f.write('VARIANT="SNPs"' + '\n')
						elif VTYPE == 'indels':
							f.write('VARIANT="Indels"' + '\n')
						f.write('sed "s/<SPECIES>/$SPECIES/g" ' + TEMPLATE_PATH + ' > ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<DATASET>|$DATASET|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<OUTDIR>|$OUTDIR|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<IN_TYPE>|$IN_TYPE|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<TE_TYPE>|$TE_TYPE|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<VTYPE>|$VTYPE|g" ' + RSCRIPT_FILE + '\n')
						f.write('sed -i "s|<VARIANT>|$VARIANT|g" ' + RSCRIPT_FILE + '\n')
						f.write('\n')
						f.write('Rscript ' + RSCRIPT_FILE + '\n')
						with open(JOBSUB_FILE, 'a+') as g:
							g.write('qsub ' + QSUB_FILE + '\n')
				else:
					print('Bad queue choice. Your only choices are hrothgar and quanah.')
