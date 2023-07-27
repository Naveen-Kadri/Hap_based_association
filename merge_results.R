##infiles  <- paste0('/cluster/work/pausch/naveen/RECOMBINATION/HAPBASED/all/male/CHR', 1:29, '/result.txt')
## chr  <- 1
infiles  <- snakemake@input[["infiles"]]
outfile  <- snakemake@output[["outfile"]]
frq_thresh  <- snakemake@wildcards [["frqthresh"]]

cs  <- rep ('numeric', 8)
cs [c(3,5)]  <- 'character'
all  <- data.frame ()
for (chr in 1:29) {
    inf  <- read.table (infiles [chr], colClasses=cs, fill=T)
    colnames (inf)  <- c('start', 'end', 'hap', 'frq', 'positions', 'myfrq',  'eff', 'se', 't', 'pval')
    inf  <- inf [!is.na (inf$pval),]
    inf  <- inf [inf$myfrq > frq_thresh,]
    ##keep only minimum at a position
    min  <- aggregate (inf$pval, by=list (inf$start),  min)
    min$chr  <- chr
    all  <- rbind (all, min)
    cat (chr, "\n ")
}


write.table (all, outfile, col.names=F,quote=F,sep="\t", row.names=F)
