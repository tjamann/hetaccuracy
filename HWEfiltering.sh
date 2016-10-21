#!bin/bash
set -e
set -u

date

KEEP=keep_seq.txt

#create the file with the hardy weinberg p values. we are going to try filtering on p values for excess hets first
vcftools --gzvcf final.GATK.break.vcf.recode.vcf.gz --hardy

#the python script to remove the variants
python filterphetexcess.py

#filter for just the inbreds and hybrids to test for heteroz. accuracy
vcftools --gzvcf allsamplesafterfiltering.mpileup.vcf.gz --keep $KEEP --remove-filtered-all --recode --out mpileup.filteredindiv
vcftools --gzvcf final.GATK.break.vcf.recode.vcf.gz --keep $KEEP --remove-filtered-all --recode --out GATK.filteredindiv

#output .012 file that serves as input for r script
bgzip mpileup.filteredindiv.recode.vcf
bgzip GATK.filteredindiv.recode.vcf
tabix -p vcf mpileup.filteredindiv.recode.vcf.gz
tabix -p vcf GATK.filteredindiv.recode.vcf.gz
vcf-merge mpileup.filteredindiv.recode.vcf.gz GATK.filteredindiv.recode.vcf.gz final.hapmap.recode.vcf.gz > mpile.GATK.true.vcf
vcftools --vcf mpile.GATK.true.vcf --012

#then run r script het.RMD
R -e "rmarkdown::render('het.Rmd')"

