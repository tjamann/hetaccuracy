#!bin/bash
set -e
set -u

date

KEEP=keep_seq.txt
GATK=final.GATK.break.vcf.recode.vcf
GATKGZ=final.GATK.break.vcf.recode.vcf.gz
PVAL=0.01
COFF=0.15

#p value is what the cutoff is for individual variants
#coff is the fraction of variants that if exceeded, the whole amplicon is removed

#to remove duplicates from bedtools join 
#awk '! a[$1" "$2]++' allsamplesafterfiltering.mpileup.vcf > allsamplesafterfiltering.dupremoved.mpileup.vcf
#the files to be compressed for the next steps

if [ -f $GATK ]
  then bgzip $GATK
  tabix -p vcf $GATKGZ
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
rm -f mpile.GATK.true.vcf
rm -f out*
rm -f mpileup*
rm -f GATK*

#filter for just the inbreds and hybrids to test for het. accuracy
vcftools --gzvcf allsamplesafterfiltering.dupremoved.mpileup.vcf.gz --keep $KEEP --remove-filtered-all --recode --out mpileup.filteredindiv
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

#mv out.012 outphetexcess.012 -i
#mv out.012.indv outphetexcess.012.indv -i
#mv out.012.pos outphetexcess.012.pos -i

cat hetaccuracies.csv 
mv hetaccuracies.csv hetaccuracies_phet.unfilt.csv
mv hetaccuraciesaverages.csv hetaccuracies_phet_averages.unfilt.csv