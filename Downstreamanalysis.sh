#!bin/bash
set -e
set -u

date

MPILEFILT=PHETE-FILTERED_allsamplesafterfiltering.dupremoved.mpileup.vcf.gz
GATKFILT=PHETE-FILTERED_final.GATK.break.vcf.recode.vcf.gz

#venn diagram
vcf-compare -a $MPILEFILT $GATKFILT final.hapmap.recode.vcf.gz > compare.txt

bgzip -d $GATKFILT
bgzip -d $MPILEFILT

#figure 4
snpEff -v -csvStats snpeffbed.txt AGPv3.27 -i bed $INTERVALS_BED -s summaryintervals > bed_targets.vcf
snpEff -v -csvStats snpeffsamtools.txt AGPv3.27 allsamplesafterfiltering.mpileup.vcf -s summarysamtools > samtools_annotation.vcf
snpEff -v -csvStats snpeffGATK.txt AGPv3.27 final.GATK.break.vcf.recode.vcf -s summaryGATK > GATK_annotation.vcf
snpEff -v -csvStats snpeffhapmap.txt AGPv3.27 final.hapmap.recode.vcf -s summaryhapmap > hapmap_annotation.vcf

#figure 5
vcftools --vcfgz $GATKFILT --plink --out sortedtasselwithnames
plink --file sortedtasselwithnames --cluster --mds-plot 2 --noweb


#FST