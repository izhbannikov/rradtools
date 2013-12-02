#rradtools.R
build_radsites <- function()#fnames,keys) 
{
   keys <- c("CCTGTG", "CTGCGA", "GAGGGA", "CATAGA", "GATCCA", "CGTTAA", "ACTGAT", "AACAAC", "GTACCG", "ATCGGG", "CATCAG", "AGTCAC", "GCACAC", "TCGTCA", "TCCACG", "TTCGAC", "ATAGTT", "GGTAAG", "GGGATT", "ACCTAA", "ACCAGT", "CTAGGC", "TCACGG", "TTCCCA")
   fnames <- c("~/bio/app/RAD/first4k.fastq")
   .Call("BuildRadSites",fnames,keys)
}