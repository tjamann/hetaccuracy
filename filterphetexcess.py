import pandas
import csv
import argparse
# SAMPLE_CALL :: python filterphetexcess.py -vcf_list allsamplesafterfiltering.mpileup.vcf final.GATK.break.vcf.recode.vcf -pval_file out.hwe --pval 0.05

def to_drop(pval_file, pval=0.05):
   tofilter = pandas.read_table(pval_file)
   loci = float(len(tofilter))
   print "NUM LOCI :: ", loci
   bonf = pval/loci
   print "CUTOFF   :: ", bonf

   drop_list = tofilter[tofilter.P_HET_EXCESS <= bonf] # Rows to call FAILHW on
   #need to filter out of vcf file
   print "TO DROP  :: ", len(drop_list)
   #save to sites to remove to dict or something and then change in vcf file those sites that fail p excess filter from PASS to FAILHWE in vcf tools
   #vcftools could then be used to remove those sites, otherwise delete this lines entirely that don't meet filter

   return drop_list[['CHR', 'POS']]

def get_lines_till_header(vcf_file):
   header_line = 0
   with open(vcf_file, 'rb') as vcffile:
       vcf_reader = csv.reader(vcffile, delimiter='\t')
       line = 0
       for row in vcf_reader:
           if len(row) != 1:
               header_line = line
               break
           line += 1

   return header_line

def update_filter_status(vcf_file, drop_list):
   # http://stackoverflow.com/questions/26896382/how-to-search-pandas-data-frame-by-index-value-and-value-in-any-column
   header_line = get_lines_till_header(vcf_file)
   vcf_pd = pandas.read_table(vcf_file, skiprows=header_line)

   # ah yes list comprehension -- the classic
   for index in drop_list.index:

       chrom = drop_list.loc[index]['CHR']
       pos = drop_list.loc[index]['POS']

       to_fail = vcf_pd[(vcf_pd['#CHROM'] == chrom ) & (vcf_pd['POS']== pos)].index
       vcf_pd.loc[to_fail, 'FILTER'] = 'FAILHW'

   # vcf_pd.set_index(['#CHROM'], inplace=True)

   return vcf_pd

def merge_old(vcf_file, updated_vcf, out_file):
   header_line = get_lines_till_header(vcf_file)
   with open(vcf_file, 'rb') as vcffile:
       first_half = list( csv.reader(vcffile, delimiter='\t') )[0:header_line]
       first_half = "\n".join([i[0] for i in first_half]) # Because the line above this returns a list of lists with only one element.
       other_half = updated_vcf.to_csv(sep='\t', index=False)

       merged_vcf = first_half + '\n' + other_half

   with open(out_file, 'wb') as merged_out:
       merged_out.write(merged_vcf)

if __name__ == '__main__':
   parser = argparse.ArgumentParser(description='Update pass/fail values on the passed .VCF file.')
   parser.add_argument('-vcf_list', nargs='*', help='List with the names of the VCF files to be processed.')
   parser.add_argument('-pval_file', help='Name of the CSV file containing calculated p-vals')
   parser.add_argument('--pval', type=float, help='Pval to use for cutoff determination')
   args = parser.parse_args()

   for vcf in args.vcf_list:
       OUT_NAME = "PHETE-FILTERED_"+ vcf
       to_drop_list = to_drop(args.pval_file, args.pval)
       updated_vcf = update_filter_status(vcf, to_drop_list)
       merge_old(vcf, updated_vcf, OUT_NAME)
