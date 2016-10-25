#!bin/bash
set -e
set -u

date

REFERENCE=/home/tjamann/Documents/LRS/extdata/Zea_mays.AGPv3.27.dna.genome.fa
GATK="java -jar /home/tjamann/bin/GenomeAnalysisTK.jar"

$GATK \
-R $REFERENCE \
-T VariantFiltration \
-o allsamplesafterfiltering.GATK.vcf \
--variant allsamples.raw.GATK.vcf \
--filterExpression "QD < 5.0" \
--filterName QDFilter \
--filterExpression "QUAL < 30.0" \
--filterName QUALFilter \
--genotypeFilterExpression "DP < 12" \
--genotypeFilterName "genotype_depth" \
--setFilteredGtToNocall

vcftools --vcf allsamplesafterfiltering.GATK.vcf --remove-filtered-all --remove-indels --maf 0.02 --recode --out final.GATK.vcf