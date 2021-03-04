import sys
import os
import argparse
import itertools
import subprocess
import fnmatch
from Bio import SeqIO

def get_args():
	parser = argparse.ArgumentParser(description="Batch script generator for MELT runs given a specific max MEI mutation rate", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--id', type=str, help='Sequence simulation ID to analyze with MELT', required=True)
	parser.add_argument('-c', '--coverage', type=int, help='Sequence simulation coverage level to analyze with MELT', required=True)
	parser.add_argument('-te', '--teseq', type=str, help='Path to TE FASTA file to analyze with MELT', required=True)
	parser.add_argument('-r', '--referencegenome', type=str, help='Path to the reference genome', required=True)
	parser.add_argument('-b', '--bam', type=str, help='Path to preprocessed BAM file to analyze with MELT', required=True)
	#required = parser.add_argument_group('required arguments')
	parser.add_argument('-m', '--rate', type=int, help="Maximum mutation rate of TE ZIP file to use", default=8)
	parser.add_argument('-np', '--proc', type=int, help='Number of cores to use for multithreaded applications', default=20)
	parser.add_argument('-od', '--outdir', type=str, help='Location of directory for output files', required=True)
	parser.add_argument('-q', '--queue', type=str, help='quanah or hrothgar', required=True)

	args = parser.parse_args()
	ID = args.id
	COV = str(args.coverage)
	TE = args.teseq
	REF = args.referencegenome
	BAM = args.bam
	RATE = args.rate
	PROC = args.proc
	OUTDIR = args.outdir
	QUEUE = args.queue
	
	return ID, COV, TE, REF, BAM, RATE, PROC, OUTDIR, QUEUE
	
ID, COV, TE, REF, BAM, RATE, PROC, OUTDIR, QUEUE = get_args()

#argument sanity checks
#if not args.telist:
#  	sys.exit('You must provide a TE list for the genome you are analyzing')

print('The simulation sequence ID is ' + ID + '.')
print('The TE sequence file is ' + TE +'.')
print('The reference "genome" is ' + REF + '.')
print('The input BAM file is ' + BAM + '.')
print('The max mutation rate for MEI ZIP file is ' + str(RATE) + ' mutations in 100 bp.')
print('Use ' + str(PROC) + ' processors.')
print('The output directory is ' + OUTDIR + '.')
print('The queue is ' + QUEUE + '.') 

BASE_DIR = "/lustre/scratch/npaulat/seq_sim"
REFDIR = "/lustre/scratch/npaulat/seq_sim/references"
#COVERAGE_DIR = os.path.join(BASE_DIR, COV)
BASE_OUTDIR = os.path.join(BASE_DIR, ID)
TE_NAME = os.path.basename(TE).split(".")[0]
TE_BED_NAME = TE_NAME + "_scaff.bed"
ZIP_NAME = TE_NAME.replace("_", "v")
ZIP_EXT = "m" + str(RATE)
ZIP_NAME = ZIP_NAME + ZIP_EXT
ZIP = ZIP_NAME + "_MELT.zip"

MUT_RATE = 'mut' + str(RATE)
MUT_RATE_DIR = os.path.join(BASE_OUTDIR, MUT_RATE)
if not os.path.exists(MUT_RATE_DIR):
	print(MUT_RATE_DIR + ' does not exist. Make it.')
	os.mkdir(MUT_RATE_DIR)
else:
	print(MUT_RATE_DIR + ' exists.')

JOBSSUBMISSON = MUT_RATE + '_' + ID + '_jobs_submission_all_' + QUEUE + '.sh'
QSUB1FILENAME = MUT_RATE + '_' + ID + '_DEL_qsub_' + QUEUE + '.sh'

DEL_DIR = 'del'
MUT_RATE_DEL_DIR = os.path.join(MUT_RATE_DIR, DEL_DIR)
if not os.path.exists(MUT_RATE_DEL_DIR):
	print(MUT_RATE_DEL_DIR + ' does not exist. Make it.')
	os.mkdir(MUT_RATE_DEL_DIR)
else:
	print(MUT_RATE_DEL_DIR + ' exists.')

#For Hrothgar queue.
if QUEUE == 'hrothgar':
	with open(QSUB1FILENAME, 'w') as f:
		f.write('#!/bin/sh' + '\n')
		f.write('#$ -V' + '\n')
		f.write('#$ -cwd' + '\n')
		f.write('#$ -S /bin/bash' + '\n')
		f.write('#$ -N ' + ID + '_DEL' + '\n')
		f.write('#$ -o $JOB_NAME.o$JOB_ID' + '\n')
		f.write('#$ -e $JOB_NAME.e$JOB_ID' + '\n')
		f.write('#$ -q Chewie' + '\n')
		f.write('#$ -pe sm ' + str(PROC) + '\n')
		f.write('#$ -P communitycluster' + '\n')
		f.write('\n')
		f.write('#The reference genome is ' + REF + '.\n')
		f.write('#The max mutation rate for MEI ZIP file is ' + str(RATE) + ' mutations in 100 bp.\n')
		f.write('#Use ' + str(PROC) + ' processors.\n')
		f.write('#The output directory is ' + MUT_RATE_DEL_DIR + '.\n')
		f.write('#The queue is ' + QUEUE + '.\n') 
		f.write('\n')
		f.write('module load intel/18.0.3.222 impi/2018.3.222 java/1.8.0 bowtie2/2.3.4 samtools/1.9\n')
		f.write('\n')
		f.write('##perl is built-in; perl5v16.3\n')
		f.write('\n')
		f.write('MELT_HOME=/lustre/work/npaulat/MELTv2.1.5\n')
		f.write('WORKDIR=' + MUT_RATE_DIR + '/del\n')
		f.write('FILEDIR=' + BASE_OUTDIR + '\n')
		f.write('REFDIR=' + REFDIR + '\n')
		f.write('\n')
		f.write('echo "Run MELT-DELETION."\n')
		f.write('\n')
		f.write('cd $WORKDIR\n')
		f.write('\n')
		f.write('echo "Begin Deletion-Genotype."\n')
		f.write('java -Xmx2G -jar $MELT_HOME/MELT.jar Deletion-Genotype -w $WORKDIR -bamfile $FILEDIR/{} -bed $REFDIR/{} -h $REFDIR/GL429795.fa\n'.format(BAM, TE_BED_NAME))
		f.write('\n')
		f.write('readlink -f $WORKDIR/*.tsv > $WORKDIR/del_list.txt\n')
		f.write('\n')
		f.write('echo "Made list of deletion.tsv (full path) files to merge into final Deletion VCF."\n')
		f.write('\n')
		f.write('echo "Begin Deletion-Merge."\n')
		f.write('java -Xmx2G -jar $MELT_HOME/MELT.jar Deletion-Merge -bed $REFDIR/{} -mergelist $WORKDIR/del_list.txt -h $REFDIR/GL429795.fa -o $WORKDIR\n'.format(TE_BED_NAME))
		f.write('\n')
		f.write('echo "' + MUT_RATE + ID + ' MELT-DELETION run completed."\n')
		with open(JOBSSUBMISSON, 'a+') as g:
			g.write('qsub ' + QSUB1FILENAME + '\n')

#For Quanah queue.
elif QUEUE == 'quanah':
	with open(QSUB1FILENAME, 'w') as f:
		f.write('#!/bin/sh' + '\n')
		f.write('#$ -V' + '\n')
		f.write('#$ -cwd' + '\n')
		f.write('#$ -S /bin/bash' + '\n')
		f.write('#$ -N ' + ID + '_DEL' + '\n')
		f.write('#$ -o $JOB_NAME.o$JOB_ID' + '\n')
		f.write('#$ -e $JOB_NAME.e$JOB_ID' + '\n')
		f.write('#$ -q omni' + '\n')
		f.write('#$ -pe sm ' + str(PROC) + '\n')
		f.write('#$ -P quanah' + '\n')
		f.write('\n')
		f.write('#The reference genome is ' + REF + '.\n')
		f.write('#The max mutation rate for MEI ZIP file is ' + str(RATE) + ' mutations in 100 bp.\n')
		f.write('#Use ' + str(PROC) + ' processors.\n')
		f.write('#The output directory is ' + MUT_RATE_DEL_DIR + '.\n')
		f.write('#The queue is ' + QUEUE + '.\n') 
		f.write('\n')
		f.write('module load intel/18.0.3.222 impi/2018.3.222 java/1.8.0 bowtie2/2.3.4 samtools/1.9\n')
		f.write('\n')
		f.write('##perl is built-in; perl5v16.3\n')
		f.write('\n')
		f.write('MELT_HOME=/lustre/work/npaulat/MELTv2.1.5\n')
		f.write('WORKDIR=' + MUT_RATE_DIR + '/del\n')
		f.write('FILEDIR=' + BASE_OUTDIR + '\n')
		f.write('REFDIR=' + REFDIR + '\n')
		f.write('\n')
		f.write('echo "Run MELT-DELETION."\n')
		f.write('\n')
		f.write('cd $WORKDIR\n')
		f.write('\n')
		f.write('echo "Begin Deletion-Genotype."\n')
		f.write('java -Xmx2G -jar $MELT_HOME/MELT.jar Deletion-Genotype -w $WORKDIR -bamfile $FILEDIR/{} -bed $REFDIR/{} -h $REFDIR/GL429795.fa\n'.format(BAM, TE_BED_NAME))
		f.write('\n')
		f.write('readlink -f $WORKDIR/*.tsv > $WORKDIR/del_list.txt\n')
		f.write('\n')
		f.write('echo "Made list of deletion.tsv (full path) files to merge into final Deletion VCF."\n')
		f.write('\n')
		f.write('echo "Begin Deletion-Merge."\n')
		f.write('java -Xmx2G -jar $MELT_HOME/MELT.jar Deletion-Merge -bed $REFDIR/{} -mergelist $WORKDIR/del_list.txt -h $REFDIR/GL429795.fa -o $WORKDIR\n'.format(TE_BED_NAME))
		f.write('\n')
		f.write('echo "' + MUT_RATE + ID + ' MELT-DELETION run completed."\n')
		with open(JOBSSUBMISSON, 'a+') as g:
			g.write('qsub ' + QSUB1FILENAME + '\n')
else:
	print('Bad queue choice. Your only choices are hrothgar and quanah.')

MUT_RATE_TE_DIR = os.path.join(MUT_RATE_DIR, ZIP_NAME)
if not os.path.exists(MUT_RATE_TE_DIR):
	print(MUT_RATE_TE_DIR + ' does not exist. Make it.')
	os.mkdir(MUT_RATE_TE_DIR)
else:
	print(MUT_RATE_TE_DIR + ' exists.')

#Create individual TE MEI ZIPs and MELT-SPLIT run qsubs, and the batch submission script
QSUB2FILENAME = MUT_RATE + '_' + ZIP_NAME + '_Single_qsub_' + QUEUE + '.sh'

#For Hrothgar queue.
if QUEUE == 'hrothgar':
	with open(QSUB2FILENAME, 'w') as h:
		h.write('#!/bin/sh' + '\n')
		h.write('#$ -V' + '\n')
		h.write('#$ -cwd' + '\n')
		h.write('#$ -S /bin/bash' + '\n')
		h.write('#$ -N ' + ID + '_Single' + '\n')
		h.write('#$ -o $JOB_NAME.o$JOB_ID' + '\n')
		h.write('#$ -e $JOB_NAME.e$JOB_ID' + '\n')
		h.write('#$ -q Chewie' + '\n')
		h.write('#$ -pe sm ' + str(PROC) + '\n')
		h.write('#$ -P communitycluster' + '\n')
		h.write('\n')
		h.write('module load intel/18.0.3.222 impi/2018.3.222 java/1.8.0 bowtie2/2.3.4 samtools/1.9\n')
		h.write('\n')
		h.write('##perl is built-in; perl5v16.3\n')
		h.write('\n')
		h.write('MELT_HOME=/lustre/work/npaulat/MELTv2.1.5\n')
		h.write('WORKDIR=' + MUT_RATE_TE_DIR + '\n')
		h.write('FILEDIR=' + BASE_OUTDIR + '\n')
		h.write('REFDIR=' + REFDIR + '\n')
		h.write('\n')
		h.write('ZIP_FILE=' + ZIP + '\n')
		h.write('\n')
		h.write('cd $WORKDIR\n')
		h.write('\n')
		h.write('#=== ' + TE + ' discovery\n')
		h.write('echo "Begin MELT-Single analysis."\n')
		h.write('java -Xmx6G -jar $MELT_HOME/MELT.jar Single -w $WORKDIR -bamfile $FILEDIR/{} -c {} -h $REFDIR/GL429795.fa -t $REFDIR/$ZIP_FILE -r 124 -n $REFDIR/scaff_genetrack.bed \n'.format(BAM, COV))
		h.write('\n')
		h.write('echo "' + TE + ' MELT-Single run completed."\n')
		with open(JOBSSUBMISSON, 'a+') as g:
			g.write('qsub ' + QSUB2FILENAME + '\n')
	
#For Quanah queue.
elif QUEUE == 'quanah':
	with open(QSUB2FILENAME, 'w') as h:
		h.write('#!/bin/sh' + '\n')
		h.write('#$ -V' + '\n')
		h.write('#$ -cwd' + '\n')
		h.write('#$ -S /bin/bash' + '\n')
		h.write('#$ -N ' + ID + '_Single' + '\n')
		h.write('#$ -o $JOB_NAME.o$JOB_ID' + '\n')
		h.write('#$ -e $JOB_NAME.e$JOB_ID' + '\n')
		h.write('#$ -q omni' + '\n')
		h.write('#$ -pe sm ' + str(PROC) + '\n')
		h.write('#$ -P quanah' + '\n')
		h.write('\n')
		h.write('module load intel/18.0.3.222 impi/2018.3.222 java/1.8.0 bowtie2/2.3.4 samtools/1.9\n')
		h.write('\n')
		h.write('##perl is built-in; perl5v16.3\n')
		h.write('\n')
		h.write('MELT_HOME=/lustre/work/npaulat/MELTv2.1.5\n')
		h.write('WORKDIR=' + MUT_RATE_TE_DIR + '\n')
		h.write('FILEDIR=' + BASE_OUTDIR + '\n')
		h.write('REFDIR=' + REFDIR + '\n')
		h.write('\n')
		h.write('ZIP_FILE=' + ZIP + '\n')
		h.write('\n')
		h.write('cd $WORKDIR\n')
		h.write('\n')
		h.write('#=== ' + TE + ' discovery\n')
		h.write('echo "Begin MELT-Single analysis."\n')
		h.write('java -Xmx6G -jar $MELT_HOME/MELT.jar Single -w $WORKDIR -bamfile $FILEDIR/{} -c {} -h $REFDIR/GL429795.fa -t $REFDIR/$ZIP_FILE -r 124 -n $REFDIR/scaff_genetrack.bed \n'.format(BAM, COV))
		h.write('\n')
		h.write('echo "' + TE + ' MELT-Single run completed."\n')
		with open(JOBSSUBMISSON, 'a+') as g:
			g.write('qsub ' + QSUB2FILENAME + '\n')
else:
	print('Bad queue choice. Your only choices are hrothgar and quanah.')
