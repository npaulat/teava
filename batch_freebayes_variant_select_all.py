#!/home/daray/anaconda3/python3

###################  INFO  ################
#  batch_pseudo.py 
#  
#  Setup script for iterative pseudoreference generation genomes or genome segments
#
#  David A Ray + Nicole S Paulat
#  v2.0 completed 8 January 2018
#  
#  This script will take in arguments from the command line and output a series of submission scripts for running on #    either Hrothgar or Quanah. The number of submission scripts will depend on the number of iterations you want 
#    to do.
#  Individual submission scripts will be deposited in a working directory, specified for each taxon or sample and  
#    the same directory as all output files from their runs.
#  A jobs submission file will be output in the main directory. Invoke this script, 'sh <jobs_submission.sh> to 
#    start the processes running. The first round will be submitted to the queue immediately. All others will wait
#    their turn.
#  
#  syntax: python batch_freebayes.py 
#           -i <# iterations> \
#           -r <path to reference genome/sample.fa> \
#           -p <prefix to identify genome/sample. will be used for output directory> \
#           -1 <path to first set of paired reads R1.fastq.gz> \
#           -2 <path to second set of paired reads R2.fastq.gz> \
#           -s <path to first set of paired reads RX_cat.fastq.gz> \ 
#           -np <# processors to use> \
#           -od <path to main working directory. all files will be output to subdirectories here> \ 
#           -q <queue to be used> \
#           -nct <# compute processors. should be < nt> \
#           -nt <# data threads. nct x nt = -np>
#           -t <variant type to call> \
#########################################################################

import sys
import os
import argparse
import itertools
import subprocess
import fnmatch
from Bio import SeqIO

def get_args():
	parser = argparse.ArgumentParser(description="Iterative pseudoreference generation with pseudo-it", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--iterations', type=int, help='Number of iterations to be performed. one iteration will not inject IUPAC ambiguities', default=7, required=True)
	parser.add_argument('-r', '--referencegenome', type=str, help='Path to the reference/contigs/scaffolds used for the first iteration', required=True)
	parser.add_argument('-p', '--prefix', help='Prefix to use on output files; this is also added to the SM field with AddOrReplaceReadGroups', required=True)
	#required = parser.add_argument_group('required arguments')
	parser.add_argument('-1', '--PE1', type=str, help="Data: PE1. PE, SE, OR PE+SE DATA IS REQUIRED", default=None, required=True)
	parser.add_argument('-2', '--PE2', type=str, help="Data: PE1. PE, SE, OR PE+SE DATA IS REQUIRED", default=None, required=True)
	parser.add_argument('-s', '--SE', type=str, help="Data: PE1. PE, SE, OR PE+SE DATA IS REQUIRED", default=None, required=True)
	parser.add_argument('-np', '--proc', type=int, help='Number of cores to use for multithreaded applications', default=1)
	parser.add_argument('-od', '--outdir', type=str, help='Location of directory for output files', required=True)
	parser.add_argument('-q', '--queue', type=str, help='quanah or hrothgar', required=True)
	parser.add_argument('-nct', '--computethreads', type=int, help="number of compute threads for the GATK's UnifiedGenotyper. total CPU usage is nct*nt", default=1)
	parser.add_argument('-nt', '--datathreads', type=int, help="number of data threads for the GATK's UnifiedGenotyper. total CPU usage is nct*nt", default=1)	
	parser.add_argument('-t', '--type', type=str, help="Variant type to call, snp or all variants. Default = snp") 

	args = parser.parse_args()
	ITER = args.iterations
	REF = args.referencegenome
	PREFIX = args.prefix
	PE1 = args.PE1
	PE2 = args.PE2
	SE = args.SE
	PROC = args.proc
	OUTDIR = args.outdir
	QUEUE = args.queue
	NCT = args.computethreads
	NT = args.datathreads
	T = args.type
	
	return PREFIX, ITER, REF, PE1, PE2, SE, PROC, OUTDIR, QUEUE, NCT, NT, T
	
PREFIX, ITER, REF, PE1, PE2, SE, PROC, OUTDIR, QUEUE, NCT, NT, T = get_args()

#argument sanity checks
#if args.iterations == 1:
#    sys.exit('You sure? Just one iteration?')
#elif not args.pe1 and not args.pe2 and not args.se:
#    sys.exit('You need to specify data with the --PE1, --PE2, and --SE options')
#elif args.pe1 and not args.pe2:
#    sys.exit('You specified PE1 but not its mate (PE2)')
#elif not args.pe1 and args.pe2:
#    sys.exit('You specified PE2 but not its mate (PE1)')
#elif not args.prefix:
#  	sys.exit('You must provide a prefix for the genome you are analyzing')

TAXONDIR = OUTDIR + '/' + PREFIX
if not os.path.exists(TAXONDIR):
	print(TAXONDIR + ' does not exist. Make it.')
	os.mkdir(TAXONDIR)
else:
	print(TAXONDIR + ' exists.')

print('Iterations = ' + str(ITER))
print('The prefix is ' + PREFIX +'.')
print('The reference genome is ' + REF + '.')
#print('The reference genome prefix is ' + REFPRE + '.')
print('The first reads are in ' + PE1 + '.')
print('The second reads are in ' + PE2 + '.')
print('The unpaired reads are in ' + SE + '.')
print('Use ' + str(PROC) + ' processors.')
print('The output directory is ' + TAXONDIR + '.')
print('The queue is ' + QUEUE + '.') 
print('The number of compute threads is ' + str(NCT) + '.')
print('The number of data threads is ' + str(NT) + '.')
if T == 'snp':
	TYPE = 'snp'
	print('Variants to be called are ' + TYPE + 's.')
elif T == 'all':
	TYPE = 'all'
	print('Variants to be called are ' + TYPE + ' variants.')
else:
	TYPE = 'snp'
	print ('No variant type indicated. Defaults to ' + TYPE + 's.')
	
#Create iteration qsubs
JOBSSUBMISSON = PREFIX + '_jobs_submission_all_' + QUEUE + '.sh'
QSUB1FILENAME = TAXONDIR + '/' + PREFIX + '_qsub_iteration1_all_' + QUEUE + '.sh'

#Iteration 1 for hrothgar.
if QUEUE == 'hrothgar':
	with open(QSUB1FILENAME, 'w') as f:
		f.write('#!/bin/sh' + '\n')
		f.write('#$ -V' + '\n')
		f.write('#$ -cwd' + '\n')
		f.write('#$ -S /bin/bash' + '\n')
		f.write('#$ -N ' + PREFIX + '.1' + '\n')
		f.write('#$ -o $JOB_NAME.o$JOB_ID' + '\n')
		f.write('#$ -e $JOB_NAME.e$JOB_ID' + '\n')
		f.write('#$ -q Chewie' + '\n')
		f.write('#$ -pe fill ' + str(PROC) + '\n')
		f.write('#$ -P communitycluster' + '\n')
		f.write('\n')
		f.write('#Iterations = ' + str(ITER) + '\n')
		f.write('#The prefix is ' + PREFIX +'.\n')
		f.write('#The reference genome is ' + REF + '.\n')
		f.write('#The first reads are in ' + PE1 + '.\n')
		f.write('#The second reads are in ' + PE2 + '.\n')
		f.write('#The unpaired reads are in ' + SE + '.\n')
		f.write('#Use ' + str(PROC) + ' processors.\n')
		f.write('#The output directory is ' + TAXONDIR + '.\n')
		f.write('#The queue is ' + QUEUE + '.\n') 
		f.write('#The number of compute threads is ' + str(NCT) + '.\n')
		f.write('#The number of data threads is ' + str(NT) + '.\n')
		f.write('#The variants to be called are ' + TYPE + '.\n')
		f.write('\n')
		f.write('module load intel bwa samtools java python')
		f.write('\n')
		f.write('cd ' + TAXONDIR + '\n')
		f.write('\n')
#		f.write('print("Create sequence dictionary for reference with picard.") \n')
		DICTPREFIX = os.path.splitext(REF)[0]
		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=' + REF + ' O=' + DICTPREFIX + '.dict TMP_DIR=tmp \n')
		f.write('\n')
#		f.write('print("Create index for reference with samtools.") \n')
		f.write('samtools faidx ' + REF + '\n')
		f.write('\n')
#		f.write('print("Create index for reference with bwa.") \n')
		f.write('bwa index ' + REF + '\n')
		f.write('\n')
#		f.write('print("Map paired-end reads with bwa mem.") \n')
		f.write('bwa mem -M -t ' + str(PROC) + ' ' + REF + ' ' + PE1 + ' ' + PE2 + ' | samtools view -Sb - > ' + PREFIX + '.iteration1.pe.bam 2> ' + PREFIX + '.iteration1.pe.bam.stderr \n')
		f.write('\n')
#		f.write('print("Map single-end reads with bwa mem.") \n')
		f.write('bwa mem -M -t ' + str(PROC) + ' ' + REF + ' ' + SE + ' | samtools view -Sb - > ' + PREFIX + '.iteration1.se.bam 2> ' + PREFIX + '.iteration1.se.bam.stderr \n')
		f.write('\n')
#		f.write('print("Merge bam files.") \n')
		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar MergeSamFiles I=' + PREFIX + '.iteration1.pe.bam I=' + PREFIX + '.iteration1.se.bam O=' + PREFIX + '.iteration1.merged.bam USE_THREADING=TRUE VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmp \n')
		f.write('\n')
#		f.write('print("Sort merged bam file.") \n')
		f.write('samtools sort -o ' + PREFIX + '.iteration1.merged.sorted.bam -T hold.sorting -@ ' + str(PROC) + ' ' + PREFIX + '.iteration1.merged.bam \n')
		f.write('\n')
#		f.write('print("Add read groups.") \n')
		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar AddOrReplaceReadGroups I=' + PREFIX +'.iteration1.merged.bam O=' + PREFIX +'.iteration1.merged.RG.bam SO=coordinate LB=' + PREFIX +'_gexome PL=illumina PU=misc SM=' + PREFIX + '  VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmp \n')
		f.write('\n')
#		f.write('print("Mark duplicates.") \n')
		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar MarkDuplicates I=' + PREFIX + '.iteration1.merged.RG.bam O=' + PREFIX + '.iteration1.merged.RG_dedup.bam VALIDATION_STRINGENCY=LENIENT M=' + PREFIX + '.iteration1.dup_metrics TMP_DIR=tmp \n')
		f.write('\n')
#		f.write('print("Index dedup.bam.") \n')
		f.write('samtools index ' + PREFIX + '.iteration1.merged.RG_dedup.bam \n')
		f.write('\n')
#		f.write('print("Create intervals list.") \n')
#		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ' + REF + ' -I ' + PREFIX + '.iteration1.merged.RG_dedup.bam -o ' + PREFIX + '.iteration1.indel_intervals.list -nt ' + str(PROC) + '\n')
#		f.write('\n')
#		f.write('print("Realign indels.") \n')
#		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T IndelRealigner -R ' + REF + ' -I ' + PREFIX + '.iteration1.merged.RG_dedup.bam -targetIntervals ' + PREFIX + '.iteration1.indel_intervals.list -o ' + PREFIX + '.iteration1.realigned.bam --filter_bases_not_stored \n')
		f.write('\n')
#		f.write('print("Create regions file.") \n')
		f.write('/lustre/work/daray/software/freebayes/scripts/fasta_generate_regions.py ' + REF + '.fai 100000 >' + REF + '.100kbp.regions.txt \n')
		f.write('\n')
#		f.write('print("Call raw variants.") \n')
		f.write('/lustre/work/daray/software/freebayes/scripts/freebayes-parallel_mod ' + REF + '.100kbp.regions.txt ' + str(PROC) + ' -f ' + REF + ' ' + PREFIX + '.iteration1.merged.RG_dedup.bam >' + PREFIX + '.iteration1.' + TYPE + '.raw.vcf \n')
		f.write('\n')
#		f.write('print("Select only SNPs.") \n')
#		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T SelectVariants -R ' + REF + ' -V ' + PREFIX + '.iteration1.raw.vcf -o ' + PREFIX + '.iteration1.snps.vcf --selectTypeToInclude SNP \n')
		f.write('\n')
#		f.write('print("Filter SNPs.") \n')
#		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T VariantFiltration -R ' + REF + ' -V ' + PREFIX + '.iteration1.snps.vcf --filterExpression "MQM < 30.0 || MQMR < 30.0 || DP < 5 || DP > 60" --filterName mqm30-mqmr30-5dp60 -o ' + PREFIX + '.iteration1.filteredsnps.vcf \n')
		f.write('\n')
#		f.write('print("Create first round snp consensus.") \n')
#		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R ' + REF + ' -o ' + PREFIX + '.iteration1.snps.consensus.fa -V ' + PREFIX + '.iteration1.filteredsnps.vcf \n')
		f.write('\n')
#		f.write('print("Filter all variants.") \n')
		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T VariantFiltration -R ' + REF + ' -V ' + PREFIX + '.iteration1.' + TYPE + '.raw.vcf --filterExpression "MQM < 30.0 || MQMR < 30.0 || DP < 5 || DP > 60" --filterName mqm30-mqmr30-5dp60 -o ' + PREFIX + '.iteration1.' + TYPE + '.filtered.vcf \n')
		f.write('\n')
#		f.write('print("Left align and trim variants.") \n')
		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T LeftAlignAndTrimVariants -R ' + REF + ' -V ' + PREFIX + '.iteration1.' + TYPE + '.filtered.vcf -o ' + PREFIX + '.iteration1.' + TYPE + '.filtered_aligned.vcf \n')
		f.write('\n')
#		f.write('print("Reduce MNPs to SNPs.") \n')
		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T VariantsToAllelicPrimitives -R ' + REF + ' -V ' + PREFIX + '.iteration1.' + TYPE + '.filtered_aligned.vcf -o ' + PREFIX + '.iteration1.' + TYPE + '.filtered_aligned_prims.vcf \n')
		f.write('\n')
#		f.write('print("Create first round all variant consensus.") \n')
		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R ' + REF + ' -o ' + PREFIX + '.iteration1.' + TYPE + '.consensus.fa -V ' + PREFIX + '.iteration1.' + TYPE + '.filtered_aligned_prims.vcf \n')
		with open(JOBSSUBMISSON, 'w') as g:
			g.write('qsub ' + QSUB1FILENAME + '\n')
	
#Other iterations for hrothgar.	
	for i in range(2,ITER + 1, 1):
		NEWREF = PREFIX + '.iteration' + str(i - 1) + '.' + TYPE + '.consensus.fa'
#		NEWREFsnp = PREFIX + '.iteration' + str(i - 1) + '.snps.consensus.fa'
		QSUBFILENAME = TAXONDIR + '/' + PREFIX + '_qsub_iteration' + str(i) + '_all_' + QUEUE + '.sh'
		with open(QSUBFILENAME, 'w') as f:
			f.write('#!/bin/sh' + '\n')
			f.write('#$ -V' + '\n')
			f.write('#$ -cwd' + '\n')
			f.write('#$ -S /bin/bash' + '\n')
			f.write('#$ -N ' + PREFIX + '.' + str(i) + '\n')
			f.write('#$ -o $JOB_NAME.o$JOB_ID' + '\n')
			f.write('#$ -e $JOB_NAME.e$JOB_ID' + '\n')
			f.write('#$ -q Chewie' + '\n')
			f.write('#$ -pe fill ' + str(PROC) + '\n')
			f.write('#$ -P communitycluster' + '\n')
			f.write('\n')
			f.write('#Iterations = ' + str(ITER) + '\n')
			f.write('#The prefix is ' + PREFIX + '.\n')
			f.write('#The reference genome is ' + REF + '.\n')
			f.write('#The first reads are in ' + PE1 + '.\n')
			f.write('#The second reads are in ' + PE2 + '.\n')
			f.write('#The unpaired reads are in ' + SE + '.\n')
			f.write('#Use ' + str(PROC) + ' processors.\n')
			f.write('#The output directory is ' + TAXONDIR + '.\n')
			f.write('#The queue is ' + QUEUE + '.\n') 
			f.write('#The number of compute threads is ' + str(NCT) + '.\n')
			f.write('#The number of data threads is ' + str(NT) + '.\n')
			f.write('#The variants to be called are ' + TYPE + '.\n')
			f.write('\n')
			f.write('module load intel bwa samtools java python')
			f.write('\n')
			f.write('cd ' + TAXONDIR + '\n')
			f.write('\n')
#			f.write('print("Create sequence dictionary for reference with picard.") \n')
			DICTPREFIX = os.path.splitext(NEWREF)[0]
			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=' + NEWREF + ' O=' + DICTPREFIX + '.dict TMP_DIR=tmp \n')
			f.write('\n')
#			f.write('print("Create index for reference with samtools.") \n')
			f.write('samtools faidx ' + NEWREF + '\n')
			f.write('\n')
#			f.write('print("Create index for reference with bwa.") \n')
			f.write('bwa index ' + NEWREF + '\n')
			f.write('\n')
#			f.write('print("Map paired-end reads with bwa mem.") \n')
			f.write('bwa mem -M -t ' + str(PROC) + ' ' + NEWREF + ' ' + PE1 + ' ' + PE2 + ' | samtools view -Sb - > ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.pe.bam 2> ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.pe.bam.stderr \n')
			f.write('\n')
#			f.write('print("Map single-end reads with bwa mem.") \n')
			f.write('bwa mem -M -t ' + str(PROC) + ' ' + NEWREF + ' ' + SE + ' | samtools view -Sb - > ' + PREFIX +'.iteration' + str(i) + '.' + TYPE + '.se.bam 2> ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.se.bam.stderr \n')
			f.write('\n')
#			f.write('print("Merge bam files.") \n')
			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar MergeSamFiles I=' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.pe.bam I=' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.se.bam O=' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.merged.bam USE_THREADING=TRUE VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmp \n')
			f.write('\n')
#			f.write('print("Sort merged bam file.") \n')
			f.write('samtools sort -o ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.merged.sorted.bam -T hold.sorting -@ ' + str(PROC) + ' ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.merged.bam \n')
			f.write('\n')
#			f.write('print("Add read groups.") \n')
			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar AddOrReplaceReadGroups I=' + PREFIX +'.iteration' + str(i) + '.' + TYPE + '.merged.bam O=' + PREFIX +'.iteration' + str(i) + '.' + TYPE + '.merged.RG.bam SO=coordinate LB=' + PREFIX + '_gexome PL=illumina PU=misc SM=' + PREFIX + '  VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmp \n')
			f.write('\n')
#			f.write('print("Mark duplicates.") \n')
			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar MarkDuplicates I=' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.merged.RG.bam O=' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.merged.RG_dedup.bam VALIDATION_STRINGENCY=LENIENT M=' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.dup_metrics TMP_DIR=tmp \n')
			f.write('\n')
#			f.write('print("Index dedup.bam.") \n')
			f.write('samtools index ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.merged.RG_dedup.bam \n')
			f.write('\n')
#			f.write('print("Create intervals list.") \n')
#			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ' + NEWREFsnp + ' -I ' + PREFIX + '.iteration' + str(i) + '.all.merged.RG_dedup.bam -o ' + PREFIX + '.iteration' + str(i) + '.all.indel_intervals.list -nt ' + str(PROC) + '\n')
#			f.write('\n')
#			f.write('print("Realign indels.") \n')
#			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T IndelRealigner -R ' + NEWREF + ' -I ' + PREFIX + '.iteration' + str(i) + '.all.merged.RG_dedup.bam -targetIntervals ' + PREFIX + '.iteration' + str(i) + '.all.indel_intervals.list -o ' + PREFIX + '.iteration' + str(i) + '.all.realigned.bam --filter_bases_not_stored \n')
#			f.write('\n')
#			f.write('print("Create regions file.") \n')
			f.write('/lustre/work/daray/software/freebayes/scripts/fasta_generate_regions.py ' + NEWREF + '.fai 100000 >' + NEWREF + '.100kbp.regions.txt \n')
			f.write('\n')
#			f.write('print("Call raw variants.") \n')
			f.write('/lustre/work/daray/software/freebayes/scripts/freebayes-parallel_mod ' + NEWREF + '.100kbp.regions.txt ' + str(PROC) + ' -f ' + NEWREF + ' ' + PREFIX + '.iteration' + str(i) + '.all.merged.RG_dedup.bam >' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.raw.vcf \n')
			f.write('\n')
#			f.write('print("Select only SNPs.") \n')
#			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T SelectVariants -R ' + NEWREFsnp + ' -V ' + PREFIX + '.iteration' + str(i) + '.snps.raw.vcf -o ' + PREFIX + '.iteration' + str(i) + '.snps.vcf --selectTypeToInclude SNP \n')
#			f.write('\n')
#			f.write('print("Filter SNPs.") \n')
#			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T VariantFiltration -R ' + NEWREFsnp + ' -V ' + PREFIX + '.iteration' + str(i) + '.snps.vcf --filterExpression "MQM < 30.0 || MQMR < 30.0 || DP < 5 || DP > 60" --filterName mqm30-mqmr30-5dp60 -o ' + PREFIX + '.iteration' + str(i) + '.filteredsnps.vcf \n')
#			f.write('\n')
#			f.write('print("Create iterative snp consensus.") \n')
#			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R ' + NEWREFsnp + ' -o ' + PREFIX + '.iteration' + str(i) + '.snps.consensus.fa -V ' + PREFIX + '.iteration' + str(i) + '.filteredsnps.vcf \n')
#			f.write('\n')
#			f.write('print("Filter all variants.") \n')
			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T VariantFiltration -R ' + NEWREF + ' -V ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.raw.vcf --filterExpression "MQM < 30.0 || MQMR < 30.0 || DP < 5 || DP > 60" --filterName mqm30-mqmr30-5dp60 -o ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.filtered.vcf \n')
			f.write('\n')
#			f.write('print("Left align and trim variants.") \n')
			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T LeftAlignAndTrimVariants -R ' + NEWREF + ' -V ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.filtered.vcf -o ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.filtered_aligned.vcf \n')
			f.write('\n')
#			f.write('print("Reduce MNPs to SNPs.") \n')
			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T VariantsToAllelicPrimitives -R ' + NEWREF + ' -V ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.filtered_aligned.vcf -o ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.filtered_aligned_prims.vcf \n')
			f.write('\n')
#			f.write('print("Create iterative all variant consensus.") \n')
			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R ' + NEWREF + ' -o ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.consensus.fa -V ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.filtered_aligned_prims.vcf \n')
			f.write('\n')
		with open(JOBSSUBMISSON, 'a') as g:
			g.write('qsub -hold_jid '+ PREFIX + '.' + str(i-1) + ' ' + QSUBFILENAME + '\n')

#Iteration 1 for quanah
elif QUEUE == 'quanah':
	with open(QSUB1FILENAME, 'w') as f:
		f.write('#!/bin/sh' + '\n')
		f.write('#$ -V' + '\n')
		f.write('#$ -cwd' + '\n')
		f.write('#$ -S /bin/bash' + '\n')
		f.write('#$ -N ' + PREFIX + '.1' + '\n')
		f.write('#$ -o $JOB_NAME.o$JOB_ID' + '\n')
		f.write('#$ -e $JOB_NAME.e$JOB_ID' + '\n')
		f.write('#$ -q omni' + '\n')
		f.write('#$ -pe sm ' + str(PROC) + '\n')
		f.write('#$ -P quanah' + '\n')
		f.write('\n')
		f.write('#Iterations = ' + str(ITER) + '\n')
		f.write('#The prefix is ' + PREFIX +'.\n')
		f.write('#The reference genome is ' + REF + '.\n')
		f.write('#The first reads are in ' + PE1 + '.\n')
		f.write('#The second reads are in ' + PE2 + '.\n')
		f.write('#The unpaired reads are in ' + SE + '.\n')
		f.write('#Use ' + str(PROC) + ' processors.\n')
		f.write('#The output directory is ' + TAXONDIR + '.\n')
		f.write('#The queue is ' + QUEUE + '.\n') 
		f.write('#The number of compute threads is ' + str(NCT) + '.\n')
		f.write('#The number of data threads is ' + str(NT) + '.\n')
		f.write('#The variants to be called are ' + TYPE + '.\n')
		f.write('\n')
		f.write('module load intel bwa samtools java python')
		f.write('\n')
		f.write('cd ' + TAXONDIR + '\n')
		f.write('\n')
#		f.write('print("Create sequence dictionary for reference with picard.") \n')
		DICTPREFIX = os.path.splitext(REF)[0]
		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=' + REF + ' O=' + DICTPREFIX + '.dict TMP_DIR=tmp \n')
		f.write('\n')
#		f.write('print("Create index for reference with samtools.") \n')
		f.write('samtools faidx ' + REF + '\n')
		f.write('\n')
#		f.write('print("Create index for reference with bwa.") \n')
		f.write('bwa index ' + REF + '\n')
		f.write('\n')
#		f.write('print("Map paired-end reads with bwa mem.") \n')
		f.write('bwa mem -M -t ' + str(PROC) + ' ' + REF + ' ' + PE1 + ' ' + PE2 + ' | samtools view -Sb - > ' + PREFIX + '.iteration1.pe.bam 2> ' + PREFIX + '.iteration1.pe.bam.stderr \n')
		f.write('\n')
#		f.write('print("Map single-end reads with bwa mem.") \n')
		f.write('bwa mem -M -t ' + str(PROC) + ' ' + REF + ' ' + SE + ' | samtools view -Sb - > ' + PREFIX + '.iteration1.se.bam 2> ' + PREFIX + '.iteration1.se.bam.stderr \n')
		f.write('\n')
#		f.write('print("Merge bam files.") \n')
		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar MergeSamFiles I=' + PREFIX + '.iteration1.pe.bam I=' + PREFIX + '.iteration1.se.bam O=' + PREFIX + '.iteration1.merged.bam USE_THREADING=TRUE VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmp \n')
		f.write('\n')
#		f.write('print("Sort merged bam file.") \n')
		f.write('samtools sort -o ' + PREFIX + '.iteration1.merged.sorted.bam -T hold.sorting -@ ' + str(PROC) + ' ' + PREFIX + '.iteration1.merged.bam \n')
		f.write('\n')
#		f.write('print("Add read groups.") \n')
		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar AddOrReplaceReadGroups I=' + PREFIX +'.iteration1.merged.bam O=' + PREFIX +'.iteration1.merged.RG.bam SO=coordinate LB=' + PREFIX +'_gexome PL=illumina PU=misc SM=' + PREFIX + '  VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmp \n')
		f.write('\n')
#		f.write('print("Mark duplicates.") \n')
		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar MarkDuplicates I=' + PREFIX + '.iteration1.merged.RG.bam O=' + PREFIX + '.iteration1.merged.RG_dedup.bam VALIDATION_STRINGENCY=LENIENT M=' + PREFIX + '.iteration1.dup_metrics TMP_DIR=tmp \n')
		f.write('\n')
#		f.write('print("Index dedup.bam.") \n')
		f.write('samtools index ' + PREFIX + '.iteration1.merged.RG_dedup.bam \n')
		f.write('\n')
#		f.write('print("Create intervals list.") \n')
#		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ' + REF + ' -I ' + PREFIX + '.iteration1.merged.RG_dedup.bam -o ' + PREFIX + '.iteration1.indel_intervals.list -nt ' + str(PROC) + '\n')
#		f.write('\n')
#		f.write('print("Realign indels.") \n')
#		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T IndelRealigner -R ' + REF + ' -I ' + PREFIX + '.iteration1.merged.RG_dedup.bam -targetIntervals ' + PREFIX + '.iteration1.indel_intervals.list -o ' + PREFIX + '.iteration1.realigned.bam --filter_bases_not_stored \n')
		f.write('\n')
#		f.write('print("Create regions file.") \n')
		f.write('/lustre/work/daray/software/freebayes/scripts/fasta_generate_regions.py ' + REF + '.fai 100000 >' + REF + '.100kbp.regions.txt \n')
		f.write('\n')
#		f.write('print("Call raw variants.") \n')
		f.write('/lustre/work/daray/software/freebayes/scripts/freebayes-parallel_mod ' + REF + '.100kbp.regions.txt ' + str(PROC) + ' -f ' + REF + ' ' + PREFIX + '.iteration1.merged.RG_dedup.bam >' + PREFIX + '.iteration1.' + TYPE + '.raw.vcf \n')
		f.write('\n')
#		f.write('print("Select only SNPs.") \n')
#		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T SelectVariants -R ' + REF + ' -V ' + PREFIX + '.iteration1.raw.vcf -o ' + PREFIX + '.iteration1.snps.vcf --selectTypeToInclude SNP \n')
		f.write('\n')
#		f.write('print("Filter SNPs.") \n')
#		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T VariantFiltration -R ' + REF + ' -V ' + PREFIX + '.iteration1.snps.vcf --filterExpression "MQM < 30.0 || MQMR < 30.0 || DP < 5 || DP > 60" --filterName mqm30-mqmr30-5dp60 -o ' + PREFIX + '.iteration1.filteredsnps.vcf \n')
		f.write('\n')
#		f.write('print("Create first round snp consensus.") \n')
#		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R ' + REF + ' -o ' + PREFIX + '.iteration1.snps.consensus.fa -V ' + PREFIX + '.iteration1.filteredsnps.vcf \n')
		f.write('\n')
#		f.write('print("Filter all variants.") \n')
		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T VariantFiltration -R ' + REF + ' -V ' + PREFIX + '.iteration1.' + TYPE + '.raw.vcf --filterExpression "MQM < 30.0 || MQMR < 30.0 || DP < 5 || DP > 60" --filterName mqm30-mqmr30-5dp60 -o ' + PREFIX + '.iteration1.' + TYPE + '.filtered.vcf \n')
		f.write('\n')
#		f.write('print("Left align and trim variants.") \n')
		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T LeftAlignAndTrimVariants -R ' + REF + ' -V ' + PREFIX + '.iteration1.' + TYPE + '.filtered.vcf -o ' + PREFIX + '.iteration1.' + TYPE + '.filtered_aligned.vcf \n')
		f.write('\n')
#		f.write('print("Reduce MNPs to SNPs.") \n')
		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T VariantsToAllelicPrimitives -R ' + REF + ' -V ' + PREFIX + '.iteration1.' + TYPE + '.filtered_aligned.vcf -o ' + PREFIX + '.iteration1.' + TYPE + '.filtered_aligned_prims.vcf \n')
		f.write('\n')
#		f.write('print("Create first round all variant consensus.") \n')
		f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R ' + REF + ' -o ' + PREFIX + '.iteration1.' + TYPE + '.consensus.fa -V ' + PREFIX + '.iteration1.' + TYPE + '.filtered_aligned_prims.vcf \n')
		with open(JOBSSUBMISSON, 'w') as g:
			g.write('qsub ' + QSUB1FILENAME + '\n')
	
#Other iterations for quanah.	
	for i in range(2,ITER + 1, 1):
		NEWREF = PREFIX + '.iteration' + str(i - 1) + '.' + TYPE + '.consensus.fa'
#		NEWREFsnp = PREFIX + '.iteration' + str(i - 1) + '.snps.consensus.fa'
		QSUBFILENAME = TAXONDIR + '/' + PREFIX + '_qsub_iteration' + str(i) + '_all_' + QUEUE + '.sh'
		with open(QSUBFILENAME, 'w') as f:
			f.write('#!/bin/sh' + '\n')
			f.write('#$ -V' + '\n')
			f.write('#$ -cwd' + '\n')
			f.write('#$ -S /bin/bash' + '\n')
			f.write('#$ -N ' + PREFIX + '.' + str(i) + '\n')
			f.write('#$ -o $JOB_NAME.o$JOB_ID' + '\n')
			f.write('#$ -e $JOB_NAME.e$JOB_ID' + '\n')
			f.write('#$ -q omni' + '\n')
			f.write('#$ -pe sm ' + str(PROC) + '\n')
			f.write('#$ -P quanah' + '\n')
			f.write('\n')
			f.write('#Iterations = ' + str(ITER) + '\n')
			f.write('#The prefix is ' + PREFIX + '.\n')
			f.write('#The reference genome is ' + REF + '.\n')
			f.write('#The first reads are in ' + PE1 + '.\n')
			f.write('#The second reads are in ' + PE2 + '.\n')
			f.write('#The unpaired reads are in ' + SE + '.\n')
			f.write('#Use ' + str(PROC) + ' processors.\n')
			f.write('#The output directory is ' + TAXONDIR + '.\n')
			f.write('#The queue is ' + QUEUE + '.\n') 
			f.write('#The number of compute threads is ' + str(NCT) + '.\n')
			f.write('#The number of data threads is ' + str(NT) + '.\n')
			f.write('#The variants to be called are ' + TYPE + '.\n')
			f.write('\n')
			f.write('module load intel bwa samtools java python')
			f.write('\n')
			f.write('cd ' + TAXONDIR + '\n')
			f.write('\n')
#			f.write('print("Create sequence dictionary for reference with picard.") \n')
			DICTPREFIX = os.path.splitext(NEWREF)[0]
			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=' + NEWREF + ' O=' + DICTPREFIX + '.dict TMP_DIR=tmp \n')
			f.write('\n')
#			f.write('print("Create index for reference with samtools.") \n')
			f.write('samtools faidx ' + NEWREF + '\n')
			f.write('\n')
#			f.write('print("Create index for reference with bwa.") \n')
			f.write('bwa index ' + NEWREF + '\n')
			f.write('\n')
#			f.write('print("Map paired-end reads with bwa mem.") \n')
			f.write('bwa mem -M -t ' + str(PROC) + ' ' + NEWREF + ' ' + PE1 + ' ' + PE2 + ' | samtools view -Sb - > ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.pe.bam 2> ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.pe.bam.stderr \n')
			f.write('\n')
#			f.write('print("Map single-end reads with bwa mem.") \n')
			f.write('bwa mem -M -t ' + str(PROC) + ' ' + NEWREF + ' ' + SE + ' | samtools view -Sb - > ' + PREFIX +'.iteration' + str(i) + '.' + TYPE + '.se.bam 2> ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.se.bam.stderr \n')
			f.write('\n')
#			f.write('print("Merge bam files.") \n')
			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar MergeSamFiles I=' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.pe.bam I=' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.se.bam O=' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.merged.bam USE_THREADING=TRUE VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmp \n')
			f.write('\n')
#			f.write('print("Sort merged bam file.") \n')
			f.write('samtools sort -o ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.merged.sorted.bam -T hold.sorting -@ ' + str(PROC) + ' ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.merged.bam \n')
			f.write('\n')
#			f.write('print("Add read groups.") \n')
			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar AddOrReplaceReadGroups I=' + PREFIX +'.iteration' + str(i) + '.' + TYPE + '.merged.bam O=' + PREFIX +'.iteration' + str(i) + '.' + TYPE + '.merged.RG.bam SO=coordinate LB=' + PREFIX + '_gexome PL=illumina PU=misc SM=' + PREFIX + '  VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmp \n')
			f.write('\n')
#			f.write('print("Mark duplicates.") \n')
			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar MarkDuplicates I=' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.merged.RG.bam O=' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.merged.RG_dedup.bam VALIDATION_STRINGENCY=LENIENT M=' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.dup_metrics TMP_DIR=tmp \n')
			f.write('\n')
#			f.write('print("Index dedup.bam.") \n')
			f.write('samtools index ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.merged.RG_dedup.bam \n')
#			f.write('\n')
#			f.write('print("Create intervals list.") \n')
#			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ' + NEWREFsnp + ' -I ' + PREFIX + '.iteration' + str(i) + '.all.merged.RG_dedup.bam -o ' + PREFIX + '.iteration' + str(i) + '.all.indel_intervals.list -nt ' + str(PROC) + '\n')
#			f.write('\n')
#			f.write('print("Realign indels.") \n')
#			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T IndelRealigner -R ' + NEWREF + ' -I ' + PREFIX + '.iteration' + str(i) + '.all.merged.RG_dedup.bam -targetIntervals ' + PREFIX + '.iteration' + str(i) + '.all.indel_intervals.list -o ' + PREFIX + '.iteration' + str(i) + '.all.realigned.bam --filter_bases_not_stored \n')
#			f.write('\n')
#			f.write('print("Create regions file.") \n')
			f.write('/lustre/work/daray/software/freebayes/scripts/fasta_generate_regions.py ' + NEWREF + '.fai 100000 >' + NEWREF + '.100kbp.regions.txt \n')
			f.write('\n')
#			f.write('print("Call raw variants.") \n')
			f.write('/lustre/work/daray/software/freebayes/scripts/freebayes-parallel_mod ' + NEWREF + '.100kbp.regions.txt ' + str(PROC) + ' -f ' + NEWREF + ' ' + PREFIX + '.iteration' + str(i) + '.all.merged.RG_dedup.bam >' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.raw.vcf \n')
#			f.write('\n')
#			f.write('print("Select only SNPs.") \n')
#			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T SelectVariants -R ' + NEWREFsnp + ' -V ' + PREFIX + '.iteration' + str(i) + '.snps.raw.vcf -o ' + PREFIX + '.iteration' + str(i) + '.snps.vcf --selectTypeToInclude SNP \n')
#			f.write('\n')
#			f.write('print("Filter SNPs.") \n')
#			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T VariantFiltration -R ' + NEWREFsnp + ' -V ' + PREFIX + '.iteration' + str(i) + '.snps.vcf --filterExpression "MQM < 30.0 || MQMR < 30.0 || DP < 5 || DP > 60" --filterName mqm30-mqmr30-5dp60 -o ' + PREFIX + '.iteration' + str(i) + '.filteredsnps.vcf \n')
#			f.write('\n')
#			f.write('print("Create iterative snp consensus.") \n')
#			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R ' + NEWREFsnp + ' -o ' + PREFIX + '.iteration' + str(i) + '.snps.consensus.fa -V ' + PREFIX + '.iteration' + str(i) + '.filteredsnps.vcf \n')
#			f.write('\n')
#			f.write('print("Filter all variants.") \n')
			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T VariantFiltration -R ' + NEWREF + ' -V ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.raw.vcf --filterExpression "MQM < 30.0 || MQMR < 30.0 || DP < 5 || DP > 60" --filterName mqm30-mqmr30-5dp60 -o ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.filtered.vcf \n')
			f.write('\n')
#			f.write('print("Left align and trim variants.") \n')
			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T LeftAlignAndTrimVariants -R ' + NEWREF + ' -V ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.filtered.vcf -o ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.filtered_aligned.vcf \n')
			f.write('\n')
#			f.write('print("Reduce MNPs to SNPs.") \n')
			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T VariantsToAllelicPrimitives -R ' + NEWREF + ' -V ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.filtered_aligned.vcf -o ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.filtered_aligned_prims.vcf \n')
			f.write('\n')
#			f.write('print("Create iterative all variant consensus.") \n')
			f.write('java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R ' + NEWREF + ' -o ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.consensus.fa -V ' + PREFIX + '.iteration' + str(i) + '.' + TYPE + '.filtered_aligned_prims.vcf \n')
			f.write('\n')
		with open(JOBSSUBMISSON, 'a') as g:
			g.write('qsub -hold_jid '+ PREFIX + '.' + str(i-1) + ' ' + QSUBFILENAME + '\n')
else:
	print('Bad queue choice. Your only choices are hrothgar and quanah.')

