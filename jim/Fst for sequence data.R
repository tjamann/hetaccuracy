#Estimate variation among counties and accessions within countries
#from sequence data at candidate genes
#data from first plate only, ARG, BOV and checks
library(hierfstat)

setwd("C:/Users/jholland/Dropbox/PLOSone_methods/Fst and clusters")
dat = read.table("VCFtotext2.txt", header = T, stringsAsFactors = T, na.strings = "N")
taxa = dat[,1]
dat = dat[,-1]

#Need to convert accession names to two heirarchical levels,
#nation and accession.
country = substr(taxa, 1, 3)
country[country %in% c("B73", "Mo1")] = "USA"
accession = sub("*_.", "", taxa)

#Make an indicator variable for which rows of data belong to checks rather than accessions
checks = ! country %in% c("ARG", "BOV")
#check what we are picking up as checks:
taxa[checks] #ok

#convert marker base pair calls to 11, 12, 22 format used by hierfstat package
for (col in 1:ncol(dat)) levels(dat[,col])[! levels(dat[,col]) %in% c("A", "C", "G", "T")] = "12"

for (col in 1:ncol(dat)) {
  if (levels(dat[,col])[1] != "12") {
    levels(dat[,col])[1] = "11"
  } else levels(dat[,col])[2] = "11"
}

for (col in 1:ncol(dat)) {
  levels(dat[,col])[!levels(dat[,col]) %in% c("11", "12")] = "22"
}

#Remove markers that are all missing in accessions (exclude last three rows)
missing_markers = apply(dat[!checks,], 2,function(x) all(is.na(x)))
any(missing_markers) #none are completely missing

# but some markers are almost all missing, like these:
summary(dat[!checks,7:8])

#let's find all markers that have > 50% missing data in the accessions and drop them
percentMissing = function(col) length(col[is.na(col)])/length(col)
percentMissing(dat[,7])

summaryMiss = apply(dat[!checks,], 2, percentMissing)
summary(summaryMiss)
keepLoci = summaryMiss < 0.5

dat2 = dat[, keepLoci]

#for hierfstat, need to convert all the variables to numeric
#for factors, just use the factor level indicator variables

country_code = as.numeric(as.factor(country))
accession_code = as.numeric(as.factor(accession))

datChr = apply(dat2, 2, as.character) #returns a character matrix
datNum = apply(datChr, c(1,2), as.numeric)#convert to numeric

dat2b = as.data.frame(cbind(country_code, accession_code, datNum))
#remove checks for Fst analysis
dat2b.nochecks = dat2b[!checks,]

#sort for hierstat package
dat2b = dat2b[order(country_code, accession_code),]


#wc() function is for weir-cockerham estimate of Fst, but it doesn't work in this case
#not sure why, something to do with missing data pattern or some high correlations among marker pairs?
#use the varcomp_global() function to get the variance components associated with each
#level of the hierarchy (countries, then accessions within countries, then plants within accessions)
fst_country_acc = varcomp.glob(levels = dat2b.nochecks[,1:2], loci = dat2b.nochecks[,3:ncol(dat2b.nochecks)], diploid = T)
fst_country_acc$overall
fst_country_acc$F

fst_country = varcomp.glob(levels = dat2b.nochecks[,1], loci = dat2b.nochecks[,3:ncol(dat2b.nochecks)], diploid = T)
fst_country$overall
fst_country$F

fst_acc = varcomp.glob(levels = dat2b.nochecks[,2], loci = dat2b.nochecks[,3:ncol(dat2b.nochecks)], diploid = T)
fst_acc$overall
fst_acc$F

#add an indicator variable for individual
#individual = seq(from = 1, to = nrow(dat2b.nochecks))
#dat2b.nochecks = cbind(individual, dat2b.nochecks)
#test.acc = test.between(data = dat2b.nochecks[,4:ncol(dat2b.nochecks)], test.lev = dat2b.nochecks$accession_code,rand.unit=dat2b.nochecks$individual,nperm=10)
#summary(test.acc)

#sort by accession code
data2b.nochecks = data2b.nochecks[order(data2b.nochecks$accession_code),]

#convert accession codes to sequential
acc_seq = as.numeric(as.factor(data2b.nochecks$accession_code))

#test.g and test.within fail on these data, I don't know why
#Gtest = test.g(data = dat2b.nochecks[,3:ncol(dat2b.nochecks)], level = acc_seq, nperm = 10)

#get the 95% bootstrap CI for F-statistics for differentiation of accessions
set.seed(1)
vc.bs = boot.vc(levels = acc_seq, loci = dat2b.nochecks[,3:ncol(dat2b.nochecks)], diploid = T, nboot = 1000 )
vc.bs$ci

#get basic summaries for accessions
overall= basic.stats(data = dat2b.nochecks[,2:ncol(dat2b.nochecks)], diploid = T, )
overall$overall

#get the weir-cockerham estimator for Fst, this is unbiased
wc(ndat = dat2b.nochecks[,2:ncol(dat2b.nochecks)], diploid = T)
