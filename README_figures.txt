used HWE_filtering_phetexcessfilter_amplicons.sh to generate vcf file.
used HWE_unfilt to compare to unfiltered dataset.
used Downstream analysis.sh to generate files for figures. 

Figure 2 representation of amplicons and ampliconrepresentation.Rmd used.the were in derived files and scripts folders respectively

Figure 3 was generated using Venn diagram 7 Nov 2016 . Rmd. 
	This requires the compare.txt file which is an output of downstreamanalysis.sh

Figure 4 was generated using the BarCharts 7 Nov 2016 . Rmd. 
	This requires compareregiontypes.csv file. 
	This was compiled using output from downstreamanalysis.sh. 
	I had to change the input csv file and the rmd script because it was misordering the legend or removing it. Most updated files are in the workdir and the fig gneration on Dropbox folder.



Figure 5. needed output from plink from donwstream analysis.sh script.
	Output from plink needed to be changed so plink.mds is now a csv file. I opened it as a txt file in excel and then saved as a csv. headings also need to be changed and merge with file that has race names (split into legends file and ids file). sorted old file with race info and new value file and then copy and pasted into old file with headers and race information. controls were updated manually with new c1 and c2 values.
	made chart using mds plot 7 nov 2016.rmd.
	
	Rmd scripts are housed in the plos one folder on Dropbox.
	
FST calculation using Jim’s FST script