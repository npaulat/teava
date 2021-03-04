#!/bin/bash

ID=$1

## Remove stats files
rm $ID.MELT.homozygote_only.vcf
rm $ID.MELT.heterozygote_only.vcf
rm $ID.MELT.homozygote_only.FP.vcf
rm $ID.MELT.homozygote_only.FN.bed
rm $ID.MELT.homozygote_only.TP.bed
rm $ID.MELT.heterozygote_only.FP.vcf
rm $ID.MELT.heterozygote_only.FN.bed
rm $ID.MELT.heterozygote_only.TP.bed
rm GT_stats.$ID.txt
rm $ID.MELT.all.TP.vcf
rm $ID.MELT.TP_no_PASS.vcf
rm TP_flag_stats.$ID.txt
rm FP_flag_stats.$ID.txt
