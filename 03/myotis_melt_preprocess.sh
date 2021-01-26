#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N myotis_preproc
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q omni
#$ -pe sm 36
#$ -P quanah

module load intel/18.0.3.222 java/1.8.0 samtools/1.9 bowtie2/2.3.4

cd /lustre/scratch/npaulat/MELTv2.1.5/combined_references/mmyo_ref_bams

LIST="mAustroriparius mBrandtii mCiliolabrum mDavidii mOccultus mSeptentrionalis_TTU mSeptentrionalis_USDA mThysanodes mVelifer mVivesi mYumanensis"

for NAME in $LIST; do samtools view -f 1 -q 20 -b -h -o ${NAME}_paired.bam ${NAME}/${NAME}.iteration1.merged.RG_dedup.bam; done

echo 'Made paired bams'

for NAME in $LIST; do samtools sort -o ${NAME}_paired.sorted.bam ${NAME}_paired.bam; samtools index ${NAME}_paired.sorted.bam; done

echo 'Sorted and indexed bams'

for i in $LIST; do java -Xmx2G -jar /lustre/work/npaulat/MELTv2.1.5/MELT.jar Preprocess -bamfile ${i}_paired.sorted.bam -h /lustre/scratch/npaulat/MELTv2.1.5/combined_references/mMyoMyo_m19_AffsNnoesSC.p1.fa; done

echo 'Preprocessed bams'
