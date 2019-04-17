#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N prep_vcfs
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q Chewie
#$ -P communitycluster
#$ -pe sm 1

module load intel vcftools bcftools/1.9

for i in Austroriparius Brandtii Ciliolabrum Davidii Occultus Sept_TTU Sept_USDA Thysanodes Velifer Vivesi Yumanensis; do vcftools --vcf m${i}.iteration1.all.raw_aligned_prims.vcf --recode --recode-INFO-all --keep-only-indels --minQ 20 --minGQ 20 --minDP 5 --max-meanDP 65 --min-alleles 2 --out ${i}_raw_align_prims_indels; bcftools filter -e 'GT="."' ${i}_raw_align_prims_indels.recode.vcf > ${i}_raw_align_prims_indels_filtGT.vcf; bcftools query -f '%CHROM\t%POS\t%END\t[%GT]\n' ${i}_raw_align_prims_indels_filtGT.vcf > ${i}_indels.vcf; done

for i in Austroriparius Brandtii Ciliolabrum Davidii Occultus Sept_TTU Sept_USDA Thysanodes Velifer Vivesi Yumanensis; do vcftools --vcf m${i}.iteration1.all.raw_aligned_prims.vcf --recode --recode-INFO-all --remove-indels --minQ 20 --minGQ 20 --minDP 5 --max-meanDP 65 --min-alleles 2 --out ${i}_raw_align_prims_snps; bcftools filter -e 'GT="."' ${i}_raw_align_prims_snps.recode.vcf > ${i}_raw_align_prims_snps_filtGT.vcf; bcftools query -f '%CHROM\t%POS\t%END\t[%GT]\n' ${i}_raw_align_prims_snps_filtGT.vcf > ${i}_snps.vcf; done

