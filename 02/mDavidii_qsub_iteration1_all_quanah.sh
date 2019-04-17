#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N mDavidii.1
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q Chewie
#$ -pe sm 20
#$ -P communitycluster

#Iterations = 1
#The prefix is mDavidii.
#The reference genome is /lustre/scratch/npaulat/freebayes_2019/myoLuc2.fa.
#The first reads are in /lustre/scratch/daray/Ray_low_cov_work/genome_data_processed/davidii/SRR534559_R1_paired.fastq.gz.
#The second reads are in /lustre/scratch/daray/Ray_low_cov_work/genome_data_processed/davidii/SRR534559_R2_paired.fastq.gz.
#The unpaired reads are in /lustre/scratch/daray/Ray_low_cov_work/genome_data_processed/davidii/SRR534559_RX_cat.fastq.gz.
#Use 36 processors.
#The output directory is /lustre/scratch/npaulat/freebayes_2019/mDavidii.
#The queue is hrothgar.
#The number of compute threads is 4.
#The number of data threads is 9.
#The variants to be called are all.

module load intel impi bwa samtools java python2
cd /lustre/scratch/npaulat/freebayes_2019/mDavidii

#java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=/lustre/scratch/npaulat/freebayes_2019/myoLuc2.fa O=/lustre/scratch/npaulat/freebayes_2019/myoLuc_MELT_scaffolds.dict TMP_DIR=tmp 

#samtools faidx /lustre/scratch/npaulat/freebayes_2019/myoLuc2.fa

#bwa index /lustre/scratch/npaulat/freebayes_2019/myoLuc2.fa

bwa mem -M -t 36 /lustre/scratch/npaulat/freebayes_2019/myoLuc2.fa /lustre/scratch/daray/Ray_low_cov_work/genome_data_processed/davidii/SRR534559_R1_paired.fastq.gz /lustre/scratch/daray/Ray_low_cov_work/genome_data_processed/davidii/SRR534559_R2_paired.fastq.gz | samtools view -Sb - > mDavidii.iteration1.pe.bam 2> mDavidii.iteration1.pe.bam.stderr 

bwa mem -M -t 36 /lustre/scratch/npaulat/freebayes_2019/myoLuc2.fa /lustre/scratch/daray/Ray_low_cov_work/genome_data_processed/davidii/SRR534559_RX_cat.fastq.gz | samtools view -Sb - > mDavidii.iteration1.se.bam 2> mDavidii.iteration1.se.bam.stderr 

java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar MergeSamFiles I=mDavidii.iteration1.pe.bam I=mDavidii.iteration1.se.bam O=mDavidii.iteration1.merged.bam USE_THREADING=TRUE VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmp 

samtools sort -o mDavidii.iteration1.merged.sorted.bam -T hold.sorting -@ 36 mDavidii.iteration1.merged.bam 

java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar AddOrReplaceReadGroups I=mDavidii.iteration1.merged.bam O=mDavidii.iteration1.merged.RG.bam SO=coordinate LB=mDavidii_gexome PL=illumina PU=misc SM=mDavidii  VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmp 

java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/picardtools/picard-tools-2.5.0/picard.jar MarkDuplicates I=mDavidii.iteration1.merged.RG.bam O=mDavidii.iteration1.merged.RG_dedup.bam VALIDATION_STRINGENCY=LENIENT M=mDavidii.iteration1.dup_metrics TMP_DIR=tmp 

samtools index mDavidii.iteration1.merged.RG_dedup.bam 

#/lustre/work/daray/software/freebayes/scripts/fasta_generate_regions.py /lustre/scratch/npaulat/freebayes_2019/myoLuc2.fa.fai 100000 >/lustre/scratch/npaulat/freebayes_2019/myoLuc2.fa.100kbp.regions.txt 

/lustre/work/daray/software/freebayes/scripts/freebayes-parallel_mod /lustre/scratch/npaulat/freebayes_2019/myoLuc2.fa.100kbp.regions.txt 36 -f /lustre/scratch/npaulat/freebayes_2019/myoLuc2.fa mDavidii.iteration1.merged.RG_dedup.bam >mDavidii.iteration1.all.raw.vcf 

java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T VariantFiltration -R /lustre/scratch/npaulat/freebayes_2019/myoLuc2.fa -V mDavidii.iteration1.all.raw.vcf --filterExpression "MQM < 30.0 || MQMR < 30.0 || DP < 5 || DP > 60" --filterName mqm30-mqmr30-5dp60 -o mDavidii.iteration1.all.filtered.vcf 

java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T LeftAlignAndTrimVariants -R /lustre/scratch/npaulat/freebayes_2019/myoLuc2.fa -V mDavidii.iteration1.all.filtered.vcf -o mDavidii.iteration1.all.filtered_aligned.vcf 

java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T VariantsToAllelicPrimitives -R /lustre/scratch/npaulat/freebayes_2019/myoLuc2.fa -V mDavidii.iteration1.all.filtered_aligned.vcf -o mDavidii.iteration1.all.filtered_aligned_prims.vcf 

java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R /lustre/scratch/npaulat/freebayes_2019/myoLuc2.fa -o mDavidii.iteration1.all.consensus.fa -V mDavidii.iteration1.all.filtered_aligned_prims.vcf 