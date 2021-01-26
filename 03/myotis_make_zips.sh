#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N MakeZIPs
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.o$JOB_ID
#$ -q Chewie
#$ -pe sm 20
#$ -P communitycluster

module load intel java samtools bowtie2

cd /lustre/scratch/npaulat/MELTv2.1.5/combined_references/zips

declare -a LIST
declare -a LIST2
LIST=(BAR1_ML HAL1-1_mMyo HAL1-1A_ML HAL1-1B_ML HAL1-1E_ML hAT1_ML hAT1_mMyo hAT2_ML hat3_ML HeliBat_N1a_ML HeliBat_N1b_ML HeliBat_N1c_ML HeliBat1 HeliBatN2 Helitron_R25_ML Helitron10_Myo Helitron12_mMyo Helitron13_mMyo Helitron-14_EF Helitron15_mMyo Helitron18_mMyo Helitron-19_EF Helitron19_mMyo Helitron2_mMyo Helitron20_mMyo Helitron22_mMyo Helitron26_mMyo Helitron-27_EF Helitron29_mMyo Helitron3_mMyo Helitron-30_EF Helitron32_mMyo Helitron-33_EF Helitron-38_EF Helitron-4_EF Helitron-43_EF Helitron5_mMyo Helitron-63_EF Helitron7_mMyo Helitron8_mMyo Helitron-9_EF L1MAB_ML L1MAB2_ML MARIN1_ML Mariner1_ML Mariner3_Ml Mariner-84_Hsal Mariner-84a_Hsal Myotis_piggyBac Myotis_Tc1_ML Myotis_Tc2_ML nhAT_186_ML nhAT17_ML nhAT1a_ML nhAT1b_ML nhAT2_730_ML nhAT2_ML nhAT-3_EF nhAT3_mMyo nhAT34_ML nhAT37_ML nhAT6_ML nhAT70_ML nMar382_Ml nMariner-5_EF nMariner-7_EF npiggy_156_ML npiggy1_mMyo npiggy111_ML npiggy165_ML npiggy2_345_ML npiggy2_41_ML npiggy2_mMyo npiggy259_ML npiggy269a_ML npiggy3_Mlyr npiggy4_mMyo piggyBac_2a_Mm piggyBac2_ML piggyBac2_Mm piggyBac2b_Mm SPIN_Ml SPIN_NA_1_Et SPIN_NA_10_Ml SPIN_NA_7_Ml SPIN_NA_8_Ml SPIN_NA_9_Ml Tc1_94_ML Tc2_122_ML Ves Ves1 Ves10 Ves11 Ves12 Ves13 Ves17 Ves2_mMyo Ves20 Ves21 Ves25 Ves27 Ves3 Ves3_ML Ves31 Ves35 Ves37 Ves38 Ves39)
LIST2=(BAR1vML HAL1v1vmMyo HAL1v1AvML HAL1v1BvML HAL1v1EvML hAT1vML hAT1vmMyo hAT2vML hat3vML HeliBatvN1avML HeliBatvN1bvML HeliBatvN1cvML HeliBat1 HeliBatN2 HelitronvR25vML Helitron10vMyo Helitron12vmMyo Helitron13vmMyo Helitronv14vEF Helitron15vmMyo Helitron18vmMyo Helitronv19vEF Helitron19vmMyo Helitron2vmMyo Helitron20vmMyo Helitron22vmMyo Helitron26vmMyo Helitronv27vEF Helitron29vmMyo Helitron3vmMyo Helitronv30vEF Helitron32vmMyo Helitronv33vEF Helitronv38vEF Helitronv4vEF Helitronv43vEF Helitron5vmMyo Helitronv63vEF Helitron7vmMyo Helitron8vmMyo Helitronv9vEF L1MABvML L1MAB2vML MARIN1vML Mariner1vML Mariner3vMl Marinerv84vHsal Marinerv84avHsal MyotisvpiggyBac MyotisvTc1vML MyotisvTc2vML nhATv186vML nhAT17vML nhAT1avML nhAT1bvML nhAT2v730vML nhAT2vML nhATv3vEF nhAT3vmMyo nhAT34vML nhAT37vML nhAT6vML nhAT70vML nMar382vMl nMarinerv5vEF nMarinerv7vEF npiggyv156vML npiggy1vmMyo npiggy111vML npiggy165vML npiggy2v345vML npiggy2v41vML npiggy2vmMyo npiggy259vML npiggy269avML npiggy3vMlyr npiggy4vmMyo piggyBacv2avMm piggyBac2vML piggyBac2vMm piggyBac2bvMm SPINvMl SPINvNAv1vEt SPINvNAv10vMl SPINvNAv7vMl SPINvNAv8vMl SPINvNAv9vMl Tc1v94vML Tc2v122vML Ves Ves1 Ves10 Ves11 Ves12 Ves13 Ves17 Ves2vmMyo Ves20 Ves21 Ves25 Ves27 Ves3 Ves3vML Ves31 Ves35 Ves37 Ves38 Ves39)

n=108

for ((i=0; i < $n; i++)); do java -Xmx1G -jar /lustre/work/npaulat/MELTv2.1.5/MELT.jar BuildTransposonZIP /lustre/scratch/npaulat/MELTv2.1.5/combined_references/fastas/${LIST[$i]}.fa /lustre/scratch/npaulat/MELTv2.1.5/combined_references/beds/${LIST[$i]}_sorted.bed ${LIST2[$i]}m10 10; done
