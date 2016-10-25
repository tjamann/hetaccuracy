#!bin/bash
set -e
set -u

date

KEEP=keep_seq.txt
PVAL=0.005

#the files to be compressed for the next steps

if [ -f final.GATK.break.vcf.recode.vcf ]
  then bgzip final.GATK.break.vcf.recode.vcf
  tabix -p vcf final.GATK.break.vcf.recode.vcf.gz
fi

if [ -f allsamplesafterfiltering.dupremoved.mpileup.vcf ]
  then bgzip allsamplesafterfiltering.dupremoved.mpileup.vcf
  tabix -p vcf allsamplesafterfiltering.dupremoved.mpileup.vcf.gz
fi

#clean up old files
rm -f PHETE-FILTERED*
rm -f mpileup.filteredindiv.recode.vcf.gz
rm -f GATK.filteredindiv.recode.vcf.gz
rm -f hetaccuracies.csv
rm -f hetaccuraciesaverages.csv
rm -f het.pdf

#create the file with the hardy weinberg p values. we are going to try filtering on p values for excess hets first
echo "GATK"
vcftools --gzvcf final.GATK.break.vcf.recode.vcf.gz --hardy  --out GATK
echo "mpileup"
vcftools --gzvcf allsamplesafterfiltering.dupremoved.mpileup.vcf.gz --hardy --out mpileup

#decompress for python script
bgzip -d final.GATK.break.vcf.recode.vcf.gz
bgzip -d allsamplesafterfiltering.dupremoved.mpileup.vcf.gz

#the python script to remove the variants. right now there is a bonferroni correction and p value of 0.05. python script should change from pass to failhw filter
python filterphetexcess.py -vcf_list allsamplesafterfiltering.dupremoved.mpileup.vcf -pval_file mpileup.hwe --pval $PVAL
python filterphetexcess.py -vcf_list final.GATK.break.vcf.recode.vcf -pval_file GATK.hwe --pval $PVAL

#recompress
#does python script generate new file names?

bgzip PHETE-FILTERED_allsamplesafterfiltering.dupremoved.mpileup.vcf
bgzip PHETE-FILTERED_final.GATK.break.vcf.recode.vcf
tabix -p vcf PHETE-FILTERED_allsamplesafterfiltering.dupremoved.mpileup.vcf.gz
tabix -p vcf PHETE-FILTERED_final.GATK.break.vcf.recode.vcf.gz

#filter for just the inbreds and hybrids to test for het. accuracy
vcftools --gzvcf PHETE-FILTERED_allsamplesafterfiltering.dupremoved.mpileup.vcf.gz --keep $KEEP --remove-filtered-all --recode --out mpileup.filteredindiv
vcftools --gzvcf PHETE-FILTERED_final.GATK.break.vcf.recode.vcf.gz --keep $KEEP --remove-filtered-all --recode --out GATK.filteredindiv

#output .012 file that serves as input for r script
bgzip mpileup.filteredindiv.recode.vcf
bgzip GATK.filteredindiv.recode.vcf
tabix -p vcf mpileup.filteredindiv.recode.vcf.gz
tabix -p vcf GATK.filteredindiv.recode.vcf.gz
vcf-merge mpileup.filteredindiv.recode.vcf.gz GATK.filteredindiv.recode.vcf.gz final.hapmap.recode.vcf.gz > mpile.GATK.true.vcf
vcftools --vcf mpile.GATK.true.vcf --012

#then run r script het.RMD
R -e "rmarkdown::render('het.Rmd')"

#mv out.012 outphetexcess.012 -i
#mv out.012.indv outphetexcess.012.indv -i
#mv out.012.pos outphetexcess.012.pos -i

cat hetaccuracies.csv 
mv hetaccuracies.csv hetaccuracies_phet020.csv