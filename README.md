How to use:

library(rradtools)

keys <- c("CCTGTG", "CTGCGA", "GAGGGA", "CATAGA", "GATCCA", "CGTTAA", "ACTGAT", "AACAAC", "GTACCG", "ATCGGG", "CATCAG", "AGTCAC", "GCACAC", "TCGTCA", "TCCACG", "TTCGAC", "ATAGTT", "GGTAAG", "GGGATT", "ACCTAA", "ACCAGT", "CTAGGC", "TCACGG", "TTCCCA")

fnames <- c("/Users/ilya/bio/app/RAD/rradtools/inst/extdata/first4k.fastq")

ans <- build_radsites(fnames, keys,astable=T,merge_sites=T,correction=F,rad_site_length=50,max_distance=3)

ans
 
