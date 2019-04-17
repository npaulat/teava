#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N <SPECIES>_GATK
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q Chewie
#$ -pe fill 20
#$ -P communitycluster

module load intel impi bwa samtools/1.3.1 java python2

cd /lustre/scratch/npaulat/freebayes_2019/<SPECIES>

java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T LeftAlignAndTrimVariants -R /lustre/scratch/npaulat/freebayes_2019/myoLuc2.fa -V <SPECIES>.iteration1.all.raw.vcf -o <SPECIES>.iteration1.all.raw_aligned.vcf

java -Xmx12g -Djava.io.tmpdir=tmp -jar /lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T VariantsToAllelicPrimitives -R /lustre/scratch/npaulat/freebayes_2019/myoLuc2.fa -V <SPECIES>.iteration1.all.raw_aligned.vcf -o <SPECIES>.iteration1.all.raw_aligned_prims.vcf 