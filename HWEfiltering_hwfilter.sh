#!bin/bash
set -e
set -u

date

KEEP=keep_seq.txt

#the files to be compressed for the next steps

if [ -f final.GATK.break.vcf.recode.vcf ]
  then bgzip final.GATK.break.vcf.recode.vcf
  tabix -p vcf final.GATK.break.vcf.recode.vcf.gz
fi

if [ -f allsamplesafterfiltering.mpileup.vcf ]
  then bgzip allsamplesafterfiltering.mpileup.vcf
  tabix -p vcf allsamplesafterfiltering.mpileup.vcf.gz
fi

#create the file with the hardy weinberg p values. we are going to try filtering on p values for excess hets first
vcftools --gzvcf final.GATK.break.vcf.recode.vcf.gz --hardy  --out GATK
vcftools --gzvcf allsamplesafterfiltering.mpileup.vcf.gz --hardy --out mpileup

#decompress for python script
bgzip -d final.GATK.break.vcf.recode.vcf.gz
bgzip -d allsamplesafterfiltering.mpileup.vcf.gz

#the python script to remove the variants. right now there is a bonferroni correction and p value of 0.05. python script should change from pass to failhw filter
python filterhw.py -vcf_list allsamplesafterfiltering.mpileup.vcf -pval_file mpileup.hwe --pval 0.01
python filterhw.py -vcf_list final.GATK.break.vcf.recode.vcf -pval_file GATK.hwe --pval 0.01

#recompress
#does python script generate new file names?

bgzip HW-FILTERED_allsamplesafterfiltering.mpileup.vcf
bgzip HW-FILTERED_final.GATK.break.vcf.recode.vcf
tabix -p vcf HW-FILTERED_allsamplesafterfiltering.mpileup.vcf.gz
tabix -p vcf HW-FILTERED_final.GATK.break.vcf.recode.vcf.gz

#filter for just the inbreds and hybrids to test for het. accuracy
vcftools --gzvcf HW-FILTERED_allsamplesafterfiltering.mpileup.vcf.gz --keep $KEEP --remove-filtered-all --recode --out mpileup.filteredindiv
vcftools --gzvcf HW-FILTERED_final.GATK.break.vcf.recode.vcf.gz --keep $KEEP --remove-filtered-all --recode --out GATK.filteredindiv

#output .012 file that serves as input for r script
bgzip mpileup.filteredindiv.recode.vcf
bgzip GATK.filteredindiv.recode.vcf
tabix -p vcf mpileup.filteredindiv.recode.vcf.gz
tabix -p vcf GATK.filteredindiv.recode.vcf.gz
vcf-merge mpileup.filteredindiv.recode.vcf.gz GATK.filteredindiv.recode.vcf.gz final.hapmap.recode.vcf.gz > HW.mpile.GATK.true.vcf
vcftools --vcf HW.mpile.GATK.true.vcf --012

#then run r script het.RMD
R -e "rmarkdown::render('het.Rmd')"

mv out.012 outhe010.012 -i
mv out.012.indv outhe010.012.indv -i
mv out.012.pos outhe010.012.pos -i
