#!/bin/bash

module load intel bedtools 

ID=$1

#create "All TP vcf file"
bedtools window -header -w 50 -a $ID.MELT.out.vcf -b $ID.insertions.sorted.bed > $ID.MELT.all.TP.vcf 

## create "All homozygote-only"-file:
grep "^#" $ID.MELT.out.vcf > $ID.MELT.homozygote_only.vcf
grep -v "^#" $ID.MELT.out.vcf | awk '$10 ~ /1\/1/' >> $ID.MELT.homozygote_only.vcf

## create "All heterozygote-only"-file:
grep "^#" $ID.MELT.out.vcf > $ID.MELT.heterozygote_only.vcf
grep -v "^#" $ID.MELT.out.vcf | awk '$10 ~ /1\/0/' >> $ID.MELT.heterozygote_only.vcf
grep -v "^#" $ID.MELT.out.vcf | awk '$10 ~ /0\/1/' >> $ID.MELT.heterozygote_only.vcf

#### Compile statistics file:
## MELT all detections (raw unfiltered output/ignoring PASS flags):
## MELT "All homozygote" (unfiltered output flagged as homozygote):
## false positives:
bedtools window -header -v -w 50 -a $ID.MELT.homozygote_only.vcf -b $ID.insertions.sorted.bed > $ID.MELT.homozygote_only.FP.vcf
#false negatives:
bedtools window -header -v -w 50 -a $ID.insertions.sorted.bed -b $ID.MELT.homozygote_only.vcf > $ID.MELT.homozygote_only.FN.bed
# true positives:
bedtools window -header -w 50 -a $ID.insertions.sorted.bed -b $ID.MELT.homozygote_only.vcf > $ID.MELT.homozygote_only.TP.bed
## MELT "All heterozygote (unfiltered insertions flagged as heterozygote):
## false positives:
bedtools window -header -v -w 50 -a $ID.MELT.heterozygote_only.vcf -b $ID.insertions.sorted.bed > $ID.MELT.heterozygote_only.FP.vcf
#false negatives:
bedtools window -header -v -w 50 -a $ID.insertions.sorted.bed -b $ID.MELT.heterozygote_only.vcf > $ID.MELT.heterozygote_only.FN.bed
# true positives:
bedtools window -header -w 50 -a $ID.insertions.sorted.bed -b $ID.MELT.heterozygote_only.vcf > $ID.MELT.heterozygote_only.TP.bed
## Generate summary TXT files:
echo -en "### MELT statistics for run $ID ####" > GT_stats.$ID.txt
echo "#### (1) All detections  ####" >> GT_stats.$ID.txt
echo "# detections:" >> GT_stats.$ID.txt
grep -v "^#" $ID.MELT.out.vcf | wc -l >> GT_stats.$ID.txt
echo "True Positives:" >> GT_stats.$ID.txt
grep -v "^#" $ID.MELT.all.TP.bed | wc -l >> GT_stats.$ID.txt
echo "False Negatives:" >> GT_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FN.bed | wc -l >> GT_stats.$ID.txt
echo "False Positives:" >> GT_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | wc -l >> GT_stats.$ID.txt
echo "#### (1a) All homozygote ####" >> GT_stats.$ID.txt
echo "# detections left:" >> GT_stats.$ID.txt
grep -v "^#" $ID.MELT.homozygote_only.vcf | wc -l >> GT_stats.$ID.txt
echo "True Positives:" >> GT_stats.$ID.txt
grep -v "^#" $ID.MELT.homozygote_only.TP.bed | wc -l >> GT_stats.$ID.txt
echo "False Positives:" >> GT_stats.$ID.txt
grep -v "^#" $ID.MELT.homozygote_only.FP.vcf | wc -l >> GT_stats.$ID.txt
echo "#### (1b) All heterozygote ####" >> GT_stats.$ID.txt
echo "# detections left:" >> GT_stats.$ID.txt
grep -v "^#" $ID.MELT.heterozygote_only.vcf | wc -l >> GT_stats.$ID.txt
echo "True Positives:" >> GT_stats.$ID.txt
grep -v "^#" $ID.MELT.heterozygote_only.TP.bed | wc -l >> GT_stats.$ID.txt
echo "False Positives:" >> GT_stats.$ID.txt
grep -v "^#" $ID.MELT.heterozygote_only.FP.vcf | wc -l >> GT_stats.$ID.txt
echo -en "\n\n\n" >> GT_stats.$ID.txt
echo "#### ESAT insertions ####" >> GT_stats.$ID.txt
grep -v "^#" $ID.insertions.sorted.bed | wc -l >> GT_stats.$ID.txt


##Compile flags stats file:
## MELT all TP detections:
grep "^#" $ID.MELT.all.TP.vcf > $ID.MELT.TP_no_PASS.vcf
grep -v "^#" $ID.MELT.all.TP.vcf | awk '$7 != "PASS"' >> $ID.MELT.TP_no_PASS.vcf

echo '#### Flags ####' >> TP_flag_stats.$ID.txt
echo "#TP w/ PASS flag" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.TP.vcf | awk '$7 == "PASS"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/o PASS flag" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.TP.vcf | awk '$7 != "PASS"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ ac0 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "ac0"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ ac0;s25 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "ac0;s25"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ hDP flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "hDP"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ hDP;ac0 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "hDP;ac0"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ hDP;ac0;s25 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "hDP;ac0;s25"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ hDP;lc flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "hDP;lc"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ hDP;lc;ac0 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "hDP;lc;ac0"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ hDP;lc;ac0;s25 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "hDP;lc;ac0;s25"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ hDP;lc;s25 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "hDP;lc;s25"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ hDP;s25 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "hDP;s25"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ lc flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "lc"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ lc;ac0 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "lc;ac0"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ lc;ac0;s25 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "lc;ac0;s25"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ lc;s25 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "lc;s25"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ rSD flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "rSD"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ rSD;ac0 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "rSD;ac0"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ rSD;ac0;s25 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "rSD;ac0;s25"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ rSD;hDP flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "rSD;hDP"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ rSD;hDP;ac0 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "rSD;hDP;ac0"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ rSD;hDP;ac0;s25 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "rSD;hDP;ac0;s25"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ rSD;hDP;lc flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "rSD;hDP;lc"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ rSD;hDP;lc;ac0 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "rSD;hDP;lc;ac0"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ rSD;hDP;lc;ac0;s25 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "rSD;hDP;lc;ac0;s25"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ rSD;hDP;lc;s25 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "rSD;hDP;lc;s25"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ rSD;hDP;s25 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "rSD;hDP;s25"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ rSD;lc flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "rSD;lc"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ rSD;lc;ac0 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "rSD;lc;ac0"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ rSD;lc;ac0;s25 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "rSD;lc;ac0;s25"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ rSD;lc;s25 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "rSD;lc;s25"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ rSD;s25 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "rSD;s25"' | wc -l >> TP_flag_stats.$ID.txt
echo "# TP w/ s25 flag only" >> TP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.TP_no_PASS.vcf | awk '$7 == "s25"' | wc -l >> TP_flag_stats.$ID.txt


## MELT all FP detections (already made)
echo '#### Flags ####' >> FP_flag_stats.$ID.txt
echo "#FP w/ PASS flag" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "PASS"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/o PASS flag" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 != "PASS"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ ac0 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "ac0"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ ac0;s25 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "ac0;s25"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ hDP flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "hDP"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ hDP;ac0 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "hDP;ac0"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ hDP;ac0;s25 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "hDP;ac0;s25"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ hDP;lc flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "hDP;lc"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ hDP;lc;ac0 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "hDP;lc;ac0"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ hDP;lc;ac0;s25 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "hDP;lc;ac0;s25"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ hDP;lc;s25 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "hDP;lc;s25"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ hDP;s25 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "hDP;s25"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ lc flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "lc"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ lc;ac0 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "lc;ac0"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ lc;ac0;s25 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "lc;ac0;s25"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ lc;s25 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "lc;s25"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ rSD flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "rSD"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ rSD;ac0 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "rSD;ac0"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ rSD;ac0;s25 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "rSD;ac0;s25"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ rSD;hDP flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "rSD;hDP"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ rSD;hDP;ac0 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "rSD;hDP;ac0"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ rSD;hDP;ac0;s25 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "rSD;hDP;ac0;s25"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ rSD;hDP;lc flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "rSD;hDP;lc"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ rSD;hDP;lc;ac0 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "rSD;hDP;lc;ac0"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ rSD;hDP;lc;ac0;s25 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "rSD;hDP;lc;ac0;s25"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ rSD;hDP;lc;s25 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "rSD;hDP;lc;s25"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ rSD;hDP;s25 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "rSD;hDP;s25"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ rSD;lc flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "rSD;lc"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ rSD;lc;ac0 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "rSD;lc;ac0"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ rSD;lc;ac0;s25 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "rSD;lc;ac0;s25"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ rSD;lc;s25 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "rSD;lc;s25"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ rSD;s25 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "rSD;s25"' | wc -l >> FP_flag_stats.$ID.txt
echo "# FP w/ s25 flag only" >> FP_flag_stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | awk '$7 == "s25"' | wc -l >> FP_flag_stats.$ID.txt
