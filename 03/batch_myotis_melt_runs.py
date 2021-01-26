import sys
import os
import argparse
import itertools
import subprocess
import fnmatch
from Bio import SeqIO

def get_args():
	parser = argparse.ArgumentParser(description="Batch script generator for MELT runs given a specific max MEI mutation rate", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-tl', '--telist', type=str, help='Path to list of TEs to analyze with MELT', required=True)
	parser.add_argument('-fl', '--falist', type=str, help='Path to list of TE FASTA file names (file basenames) to analyze with MELT, must be in same order as TE list', required=True)
	parser.add_argument('-zl', '--ziplist', type=str, help='Path to list of TE ZIP file names (basenames) to analyze with MELT, must be in same order as TE list', required=True)
	parser.add_argument('-r', '--referencegenome', type=str, help='Path to the reference genome', required=True)
	#required = parser.add_argument_group('required arguments')
	parser.add_argument('-m', '--rate', type=int, help="Maximum mutation rate for each TE MEI ZIP file to be made", default=5)
	parser.add_argument('-np', '--proc', type=int, help='Number of cores to use for multithreaded applications', default=12)
	parser.add_argument('-od', '--outdir', type=str, help='Location of directory for output files', required=True)
	parser.add_argument('-q', '--queue', type=str, help='quanah or hrothgar', required=True)

	args = parser.parse_args()
	TE_LIST = args.telist
	FASTA_LIST = args.falist
	ZIP_LIST = args.ziplist
	REF = args.referencegenome
	RATE = args.rate
	PROC = args.proc
	OUTDIR = args.outdir
	QUEUE = args.queue
	
	return TE_LIST, FASTA_LIST, ZIP_LIST, REF, RATE, PROC, OUTDIR, QUEUE
	
TE_LIST, FASTA_LIST, ZIP_LIST, REF, RATE, PROC, OUTDIR, QUEUE = get_args()

#argument sanity checks
#if not args.telist:
#  	sys.exit('You must provide a TE list for the genome you are analyzing')

print('The TE list is ' + TE_LIST +'.')
print('The TE FASTA file list is ' + FASTA_LIST +'.')
print('The TE ZIP file list is ' + ZIP_LIST +'.')
print('The reference genome is ' + REF + '.')
print('The max mutation rate for MEI ZIP file is ' + str(RATE) + ' mutations in 100 bp.')
print('Use ' + str(PROC) + ' processors.')
print('The output directory is ' + OUTDIR + '.')
print('The queue is ' + QUEUE + '.') 

MUT_RATE = 'mut' + str(RATE)
MUT_RATE_DIR = os.path.join(OUTDIR, MUT_RATE)
if not os.path.exists(MUT_RATE_DIR):
	print(MUT_RATE_DIR + ' does not exist. Make it.')
	os.mkdir(MUT_RATE_DIR)
else:
	print(MUT_RATE_DIR + ' exists.')

JOBSSUBMISSON = MUT_RATE + '_jobs_submission_all_' + QUEUE + '.sh'
QSUB1FILENAME = MUT_RATE + '_DEL_qsub_' + QUEUE + '.sh'

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
		f.write('#$ -N ' + MUT_RATE + '_DEL' + '\n')
		f.write('#$ -o $JOB_NAME.o$JOB_ID' + '\n')
		f.write('#$ -e $JOB_NAME.e$JOB_ID' + '\n')
		f.write('#$ -q Chewie' + '\n')
		f.write('#$ -pe sm ' + str(PROC) + '\n')
		f.write('#$ -P communitycluster' + '\n')
		f.write('\n')
		f.write('#The reference genome is ' + REF + '.\n')
		f.write('#The max mutation rate for MEI ZIP file is ' + str(RATE) + ' mutations in 100 bp.\n')
		f.write('#Use ' + str(PROC) + ' processors.\n')
		f.write('#The output directory is ' + MUT_RATE_DIR + '/del.\n')
		f.write('#The queue is ' + QUEUE + '.\n') 
		f.write('\n')
		f.write('module load intel/18.0.3.222 impi/2018.3.222 java/1.8.0 bowtie2/2.3.4 samtools/1.9\n')
		f.write('\n')
		f.write('##perl is built-in; perl5v16.3\n')
		f.write('\n')
		f.write('MELT_HOME=/lustre/work/npaulat/MELTv2.1.5\n')
		f.write('WORKDIR=' + MUT_RATE_DIR + '/del\n')
		#f.write('FILEDIR=/lustre/scratch/npaulat/MELTv2.1.5/combined_references\n')
		f.write('FILEDIR=/lustre/scratch/npaulat/MELT/combined_references\n')
		f.write('\n')
		f.write('echo "Run MELT-DELETION."\n')
		f.write('\n')
		f.write('cd $WORKDIR\n')
		f.write('\n')
		f.write('echo "Begin Deletion-Genotype."\n')
		f.write('for i in mAustroriparius mBrandtii mCiliolabrum mDavidii mOccultus mSeptentrionalis_TTU mSeptentrionalis_USDA mThysanodes mVelifer mVivesi mYumanensis; do java -Xmx2G -jar $MELT_HOME/MELT.jar Deletion-Genotype -w $WORKDIR -bamfile $FILEDIR/bams/${i}_paired.sorted.bam -bed $FILEDIR/beds/all_TEs_filtered.bed -h ' + REF + '; done\n')
		f.write('\n')
		f.write('readlink -f $WORKDIR/*.tsv > $WORKDIR/del_list.txt\n')
		f.write('\n')
		f.write('echo "Made list of deletion.tsv (full path) files to merge into final Deletion VCF."\n')
		f.write('\n')
		f.write('echo "Begin Deletion-Merge."\n')
		f.write('java -Xmx2G -jar $MELT_HOME/MELT.jar Deletion-Merge -bed $FILEDIR/beds/all_TEs_filtered.bed -mergelist $WORKDIR/del_list.txt -h ' + REF + ' -o $WORKDIR\n')
		f.write('\n')
		f.write('echo "' + MUT_RATE + ' MELT-DELETION run completed."\n')
		with open(JOBSSUBMISSON, 'a+') as g:
			g.write('qsub ' + QSUB1FILENAME + '\n')

#For Quanah queue.
elif QUEUE == 'quanah':
	with open(QSUB1FILENAME, 'w') as f:
		f.write('#!/bin/sh' + '\n')
		f.write('#$ -V' + '\n')
		f.write('#$ -cwd' + '\n')
		f.write('#$ -S /bin/bash' + '\n')
		f.write('#$ -N ' + MUT_RATE + '_DEL' + '\n')
		f.write('#$ -o $JOB_NAME.o$JOB_ID' + '\n')
		f.write('#$ -e $JOB_NAME.e$JOB_ID' + '\n')
		f.write('#$ -q omni' + '\n')
		f.write('#$ -pe sm ' + str(PROC) + '\n')
		f.write('#$ -P quanah' + '\n')
		f.write('\n')
		f.write('#The reference genome is ' + REF + '.\n')
		f.write('#The max mutation rate for MEI ZIP file is ' + str(RATE) + ' mutations in 100 bp.\n')
		f.write('#Use ' + str(PROC) + ' processors.\n')
		f.write('#The output directory is ' + MUT_RATE_DIR + '/del.\n')
		f.write('#The queue is ' + QUEUE + '.\n') 
		f.write('\n')
		f.write('module load intel/18.0.3.222 impi/2018.3.222 java/1.8.0 bowtie2/2.3.4 samtools/1.9\n')
		f.write('\n')
		f.write('##perl is built-in; perl5v16.3\n')
		f.write('\n')
		f.write('MELT_HOME=/lustre/work/npaulat/MELTv2.1.5\n')
		f.write('WORKDIR=' + MUT_RATE_DIR + '/del\n')
		#f.write('FILEDIR=/lustre/scratch/npaulat/MELTv2.1.5/combined_references\n')
		f.write('FILEDIR=/lustre/scratch/npaulat/MELT/combined_references\n')
		f.write('\n')
		f.write('echo "Run MELT-DELETION."\n')
		f.write('\n')
		f.write('cd $WORKDIR\n')
		f.write('\n')
		f.write('echo "Begin Deletion-Genotype."\n')
		f.write('for i in mAustroriparius mBrandtii mCiliolabrum mDavidii mOccultus mSeptentrionalis_TTU mSeptentrionalis_USDA mThysanodes mVelifer mVivesi mYumanensis; do java -Xmx2G -jar $MELT_HOME/MELT.jar Deletion-Genotype -w $WORKDIR -bamfile $FILEDIR/bams/${i}_paired.sorted.bam -bed $FILEDIR/beds/all_TEs_filtered.bed -h ' + REF + '; done\n')
		f.write('\n')
		f.write('readlink -f $WORKDIR/*.tsv > $WORKDIR/del_list.txt\n')
		f.write('\n')
		f.write('echo "Made list of deletion.tsv (full path) files to merge into final Deletion VCF."\n')
		f.write('\n')
		f.write('echo "Begin Deletion-Merge."\n')
		f.write('java -Xmx2G -jar $MELT_HOME/MELT.jar Deletion-Merge -bed $FILEDIR/beds/all_TEs_filtered.bed -mergelist $WORKDIR/del_list.txt -h ' + REF + ' -o $WORKDIR\n')
		f.write('\n')
		f.write('echo "' + MUT_RATE + ' MELT-DELETION run completed."\n')
		with open(JOBSSUBMISSON, 'a+') as g:
			g.write('qsub ' + QSUB1FILENAME + '\n')
else:
	print('Bad queue choice. Your only choices are hrothgar and quanah.')

#with open("/lustre/scratch/npaulat/MELTv2.1.5/references/te_list_may.txt", "r") as d:
with open(TE_LIST, "r") as d:
	TES = d.read().split(" ")

#with open("/lustre/scratch/npaulat/MELTv2.1.5/references/te_fasta_names_may.txt", "r") as d:
with open(FASTA_LIST, "r") as d:
	FASTAS = d.read().split(" ")

#with open("/lustre/scratch/npaulat/MELTv2.1.5/references/zip_te_names_may.txt", "r") as d:
with open(ZIP_LIST, "r") as d:
	ZIPS = d.read().split(" ")

#for TE in TES:
for TE, FASTA, ZIP in zip(TES, FASTAS, ZIPS):
	MUT_RATE_TE_DIR = os.path.join(MUT_RATE_DIR, ZIP)
	if not os.path.exists(MUT_RATE_TE_DIR):
		print(MUT_RATE_TE_DIR + ' does not exist. Make it.')
		os.mkdir(MUT_RATE_TE_DIR)
	else:
		print(MUT_RATE_TE_DIR + ' exists.')

	#Create individual TE MEI ZIPs and MELT-SPLIT run qsubs, and the batch submission script
	QSUB2FILENAME = MUT_RATE + '_' + ZIP + '_SPLIT_qsub_' + QUEUE + '.sh'

	#For Hrothgar queue.
	if QUEUE == 'hrothgar':
		with open(QSUB2FILENAME, 'w') as h:
			h.write('#!/bin/sh' + '\n')
			h.write('#$ -V' + '\n')
			h.write('#$ -cwd' + '\n')
			h.write('#$ -S /bin/bash' + '\n')
			h.write('#$ -N ' + ZIP + '_' + str(RATE) + '_SPLIT' + '\n')
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
			h.write('WORKDIR=' + MUT_RATE_DIR + '\n')
			#h.write('FILEDIR=/lustre/scratch/npaulat/MELTv2.1.5/combined_references\n')
			h.write('FILEDIR=/lustre/scratch/npaulat/MELT/combined_references\n')
#			h.write('\n')
#			h.write('cd $FILEDIR/zips\n')
			h.write('\n')
			h.write('ZIP_NAME=' + ZIP + 'm' + str(RATE) + '\n')
			h.write('ZIP_FILE=$ZIP_NAME"_MELT.zip"\n')
			h.write('\n')
#			h.write('echo "Create ' + ZIP + ' MEI ZIP with mutation rate max of ' + str(RATE) + ' reference file."\n')
#			h.write('java -Xmx1G -jar $MELT_HOME/MELT.jar BuildTransposonZIP $FILEDIR/fastas/' + FASTA + '.fa' + ' $FILEDIR/beds/' + TE + '.bed $ZIP_NAME ' + str(RATE) + '\n')
#			h.write('\n')
			h.write('cd $WORKDIR\n')
			h.write('\n')
			h.write('#=== ' + TE + ' discovery\n')
			h.write('echo "Begin IndivAnalysis."\n')
			h.write('for i in mAustroriparius mBrandtii mCiliolabrum mDavidii mOccultus mSeptentrionalis_TTU mSeptentrionalis_USDA mThysanodes mVelifer mVivesi mYumanensis; do java -Xmx6G -jar $MELT_HOME/MELT.jar IndivAnalysis -w $WORKDIR/' + ZIP + ' -bamfile $FILEDIR/bams/${i}_paired.sorted.bam -c 14 -h ' + REF + ' -t $FILEDIR/zips/$ZIP_FILE -r 150; done\n')
			h.write('\n')
			h.write('echo "Begin GroupAnalysis."\n')
			h.write('java -Xmx6G -jar $MELT_HOME/MELT.jar GroupAnalysis -discoverydir $WORKDIR/' + ZIP + ' -h ' + REF + ' -n $FILEDIR/mMyo_empty_annot.bed -t $FILEDIR/zips/$ZIP_FILE -w $WORKDIR/' + ZIP + ' -r 150\n')
			h.write('\n')
			h.write('echo "Begin Genotyping."\n')
			h.write('for i in mAustroriparius mBrandtii mCiliolabrum mDavidii mOccultus mSeptentrionalis_TTU mSeptentrionalis_USDA mThysanodes mVelifer mVivesi mYumanensis; do java -Xmx6G -jar $MELT_HOME/MELT.jar Genotype -w $WORKDIR/' + ZIP + ' -bamfile $FILEDIR/bams/${i}_paired.sorted.bam -h ' + REF + ' -t $FILEDIR/zips/$ZIP_FILE -p $WORKDIR/' + ZIP + '; done\n')
			h.write('\n')
			h.write('echo "Generate mei list from .tsv files."\n')
			h.write('ls $WORKDIR/' + ZIP + '/*.tsv > $WORKDIR/' + ZIP + '/mei_list.txt\n')
			h.write('\n')
			h.write('echo "Begin MakeVCF."\n')
			h.write('java -Xmx6G -jar $MELT_HOME/MELT.jar MakeVCF -genotypingdir $WORKDIR/' + ZIP + ' -h ' + REF + ' -t $FILEDIR/zips/$ZIP_FILE -w $WORKDIR/' + ZIP + ' -p $WORKDIR/' + ZIP + '\n')
			h.write('\n')
			h.write('echo "' + TE + ' MELT-SPLIT run completed."\n')
			with open(JOBSSUBMISSON, 'a+') as g:
				g.write('qsub ' + QSUB2FILENAME + '\n')
	
	#For Quanah queue.
	elif QUEUE == 'quanah':
		with open(QSUB2FILENAME, 'w') as h:
			h.write('#!/bin/sh' + '\n')
			h.write('#$ -V' + '\n')
			h.write('#$ -cwd' + '\n')
			h.write('#$ -S /bin/bash' + '\n')
			h.write('#$ -N ' + TE + '_' + str(RATE) + '_SPLIT' + '\n')
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
			h.write('WORKDIR=' + MUT_RATE_DIR + '\n')
			#h.write('FILEDIR=/lustre/scratch/npaulat/MELTv2.1.5/combined_references\n')
			h.write('FILEDIR=/lustre/scratch/npaulat/MELT/combined_references\n')
#			h.write('\n')
#			h.write('cd $FILEDIR/zips\n')
			h.write('\n')
			h.write('ZIP_NAME=' + ZIP + 'm' + str(RATE) + '\n')
			h.write('ZIP_FILE=$ZIP_NAME"_MELT.zip"\n')
			h.write('\n')
#			h.write('echo "Create ' + ZIP + ' MEI ZIP with mutation rate max of ' + str(RATE) + ' reference file."\n')
#			h.write('java -Xmx1G -jar $MELT_HOME/MELT.jar BuildTransposonZIP $FILEDIR/fastas/' + FASTA + '.fa' + ' $FILEDIR/beds/' + TE + '.bed $ZIP_NAME ' + str(RATE) + '\n')
#			h.write('\n')
			h.write('cd $WORKDIR\n')
			h.write('\n')
			h.write('#=== ' + TE + ' discovery\n')
			h.write('echo "Begin IndivAnalysis."\n')
			h.write('for i in mAustroriparius mBrandtii mCiliolabrum mDavidii mOccultus mSeptentrionalis_TTU mSeptentrionalis_USDA mThysanodes mVelifer mVivesi mYumanensis; do java -Xmx6G -jar $MELT_HOME/MELT.jar IndivAnalysis -w $WORKDIR/' + ZIP + ' -bamfile $FILEDIR/bams/${i}_paired.sorted.bam -c 14 -h ' + REF + ' -t $FILEDIR/zips/$ZIP_FILE -r 150; done\n')
			h.write('\n')
			h.write('echo "Begin GroupAnalysis."\n')
			h.write('java -Xmx6G -jar $MELT_HOME/MELT.jar GroupAnalysis -discoverydir $WORKDIR/' + ZIP + ' -h $FILEDIR/myoLuc2.fa -n $FILEDIR/mMyo_empty_annot.bed -t $FILEDIR/zips/$ZIP_FILE -w $WORKDIR/' + ZIP + ' -r 150\n')
			h.write('\n')
			h.write('echo "Begin Genotyping."\n')
			h.write('for i in mAustroriparius mBrandtii mCiliolabrum mDavidii mOccultus mSeptentrionalis_TTU mSeptentrionalis_USDA mThysanodes mVelifer mVivesi mYumanensis; do java -Xmx6G -jar $MELT_HOME/MELT.jar Genotype -w $WORKDIR/' + ZIP + ' -bamfile $FILEDIR/bams/${i}_paired.sorted.bam -h ' + REF + ' -t $FILEDIR/zips/$ZIP_FILE -p $WORKDIR/' + ZIP + '; done\n')
			h.write('\n')
			h.write('echo "Generate mei list from .tsv files."\n')
			h.write('ls $WORKDIR/' + ZIP + '/*.tsv > $WORKDIR/' + ZIP + '/mei_list.txt\n')
			h.write('\n')
			h.write('echo "Begin MakeVCF."\n')
			h.write('java -Xmx6G -jar $MELT_HOME/MELT.jar MakeVCF -genotypingdir $WORKDIR/' + ZIP + ' -h ' + REF + ' -t $FILEDIR/zips/$ZIP_FILE -w $WORKDIR/' + ZIP + ' -p $WORKDIR/' + ZIP + '\n')
			h.write('\n')
			h.write('echo "' + TE + ' MELT-SPLIT run completed."\n')
			with open(JOBSSUBMISSON, 'a+') as g:
				g.write('qsub ' + QSUB2FILENAME + '\n')
	else:
		print('Bad queue choice. Your only choices are hrothgar and quanah.')
