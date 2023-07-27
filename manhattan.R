#rm (list = ls (all=T))
##infile  <-  '/cluster/work/pausch/naveen/RECOMBINATION/GWAS/all/male/result.mlma'

infile  <- snakemake@input[["infile"]]
plotfile  <- snakemake@output[["plotfile"]]

#-------------------------------------------------
res <- read.table (infile, head=T)
chrcol <- 3
poscol <- 1
pcol <- 2
#frqcol  <- 6
colors <- c("darkblue", "slategray")
fontsize=1.5  #font size for the ylab
#leave larger margins for ylab
mymar <- c (4, 6, 2, 2)

pdf (plotfile, height=12,width=20)
par (mar=mymar, cex.lab=1.5,cex.axis=1.5)
thresh <- 6
#-------------------------------------------------

#alternate the color for consecutive chromosomes
mycol = rep (NA, nrow (res))
mycol [res[, chrcol] %%2 == 0] = colors [1]
mycol [res[, chrcol] %%2 == 1] = colors [2]


#get the chromosome sizes
chrs <- sort(unique (res [, chrcol]))
sizes <- c ()
for (chr in 1:length (chrs)) {
    sizes [chr] <- max(res [res [, chrcol] == chr , poscol])
    cat (chr ,  "\n")
}


#continuous positions
toadd <- c (0, cumsum (as.numeric (sizes )))
toadd <- toadd [-c(length(toadd))]

res$pos <- res[, poscol] + toadd [res [, chrcol]]


#res  <- res [res[,frqcol] > 0 ,]
maxi <- max (-log10(res[, pcol])) ; maxi
plot (res$pos, -log10(res[, pcol]), col=mycol, xaxt="n", ylab=expression (-log [10](italic("P"))), ylim =c(0, maxi*1.1), cex.main=1 , cex.lab=fontsize, cex.axis=1, xlab="", pch=20, cex=1.5,yaxt="n")
axis (side=2, las=2)

#position for chr label.. midpoint
ats <- (sizes/2) + toadd

#axis usually no space for axis
#axis (side=1,  at=ats, labels=1:29)
text (ats, rep(maxi*1.1, length(ats)), labels=chrs, cex=1)

#vertical lines separting the chromosomes
abline (v= sizes + toadd,col='gray', lty=2 )

#threshold for significance
abline (h =thresh,  col='red', lty=2 ,lwd=2)

dev.off ()
