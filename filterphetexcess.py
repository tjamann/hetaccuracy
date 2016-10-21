import pandas
tofilter = pandas.read_table('out.hwe')
loci = float(len(tofilter))
print(loci)
pval = float(0.05)
bonf = pval/loci
print(bonf)
filteredoutexcess = tofilter[tofilter.P_HET_EXCESS >= bonf] # these are the sites i want to keep
len(filteredoutexcess) #to show how many lines/variants new file should have
#need to filter out of vcf file
#save to sites to remove to dict or something and then change in vcf file those sites that fail p excess filter from PASS to FAILHWE in vcf tools
#vcftools could then be used to remove those sites, otherwise delete this lines entirely that don't meet filter