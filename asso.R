args <- commandArgs(trailingOnly=T)
startpos=args [1]
end=args [2]
hap=args [3]
freq =args [4]
positions  <- args [5]
info <- args [1:5]
resfile  <- args [6]
tmpfile  <- args [7]
npc <- as.numeric (args [8])
cs  <- rep ('numeric', npc+3)
cs [1]  <- 'character'
inf <- read.table (tmpfile,sep="\t", colClasses=cs)


colnames (inf)  <- c("id", "hap", paste0('pc', 1:npc),"pheno")

inf <- inf [ , -c(1)]
mymod <- lm (inf$pheno ~ ., data=inf)
summary  <- coef (summary (mymod))
mycol  <- which (row.names (summary) =="hap")
myres  <- as.numeric (summary [mycol,])
if (length (myres) == 0) {
    myres  <- rep(NA, 4 )
}

frq  <- mean (inf$hap)/2
out  <- c (info,frq,myres)
write(out, resfile, append=T, ncolumns=10)

