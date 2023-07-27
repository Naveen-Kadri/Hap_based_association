grr_file  <- snakemake@input [["grr_file"]]
good_gams  <- snakemake@input [["good_gams"]]
sex  <- snakemake@wildcards [["sex"]]
outfile_raw  <- snakemake@output[["outfile_raw"]]
outfile_mean  <- snakemake@output[["outfile_mean"]]
breed <- snakemake@wildcards[["breed"]]

#grr_file  <- '/cluster/work/pausch/naveen/RECOMBINATION/REFALT/LINKPHASE_RUN2/bv/cleaned_GRR.txt'
#good_gams  <- '/cluster/work/pausch/naveen/RECOMBINATION/REFALT/LINKPHASE_RUN2/bv/good_gams.txt'
#sex  <- 'female'

if (sex == "female") {
    sex  <- 2
}else {
    sex  <- 1
}


good_gams  <- read.table (good_gams, head=T)
good_gams  <-   unique (good_gams [good_gams$sex == sex , 1])

grr  <- read.table (grr_file, head=T)
grr  <- grr [grr$gam %in%  good_gams, c('par', 'GRR')]

##remove outliers
mymean  <- mean (grr$GRR)
mysd  <- sd (grr$GRR)
std  <-  abs ( (grr$GRR - mymean ) / mysd)
grr  <- grr [std <=5,]

out_raw = grr [, c(1,1,2)]


means  <- aggregate (grr$GRR, by=list (grr$par),mean)
out_mean  <- means [,c(1,1,2)]



##add breed info
out_mean$breed  <- breed
out_raw$breed  <- breed

write.table (out_raw, outfile_raw, col.names=F, row.names=F,quote=F,sep="\t")
write.table (out_mean, outfile_mean, col.names=F, row.names=F,quote=F,sep="\t")
