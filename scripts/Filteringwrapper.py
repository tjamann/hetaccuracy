import argparse
import os

# Default values
HWEFILTER_NAME = "sh HWEfiltering_phetexcessfilter.sh"
KEEP = "keep_seq.txt"
GATK = "final.GATK.vcf.recode.vcf"
GATKGZ = "final.GATK.vcf.recode.vcf.gz"
PVAL = 0.001
COFF = 0.25


parser = argparse.ArgumentParser(description='Calls HWEfilter with the following arguments')
parser.add_argument('-KEEP', const=KEEP, default=KEEP, nargs='?', help='Name of the keep file')
parser.add_argument('-GATK', const=GATK, default=GATK, nargs='?', help='Name of the GATK file')
parser.add_argument('-GATKGZ', const=GATKGZ, default=GATKGZ, nargs='?', help='Name of the GATKGZ file')
parser.add_argument('-PVAL', const=PVAL, default=PVAL, nargs='?', type=float, help='Pval to use for cutoff determination')
parser.add_argument('-COFF', const=COFF, default=COFF, nargs='?', type=float, help='% of bad vairatiants on an amplicon, if greater than this then throws all variants on amplicon out')

args = parser.parse_args()

call_str = "{} {} {} {} {} {}".format(
    HWEFILTER_NAME,
    args.KEEP,
    args.GATK,
    args.GATKGZ,
    args.PVAL,
    args.COFF
)


os.system(call_str)