##infile  <- '/cluster/work/pausch/naveen/RECOMBINATION/REFALT/GRM/male/genome.eigenvec'

infile  <- snakemake@input[["infile"]]
plotfile  <- snakemake@output[["plotfile"]]

inf  <- read.table (infile, head=F)
inf$breed  <- substr (inf [,2],1,2)
breeds  <- c('bv','fv')
cols  <- c('violet', 'orange')
inf$col  <- cols [1]
inf [inf$breed == breeds [2], "col" ]  <- cols [2] 
pdf (plotfile,height=6, width=8)
par (mar = c (6,6,4,4),cex.lab=1.5, cex.axis=1.5, las=1,mgp=c(4,1,0))
xlabs  <- 'PC1'
ylabs  <- 'PC2'
plot (inf [,3:4],  xlab=xlabs, ylab=ylabs, pch=21, col='gray', bg=inf$col,cex=1.5)
ns  <-  c(sum(inf$breed == breeds [1]), sum(inf$breed == breeds [2]))
ltext  <- paste0(breeds, " ( ", ns, " ) "); ltext

legend ('bottom', horiz=T,legend=ltext, col=cols, bty='n', cex=1.8, pch=19)
dev.off ()
