#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file)", call.=FALSE)
} else if (length(args)==1) {
  # default output name
  args[2] = "kmers.pdf"
}

if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2",repos='http://cran.us.r-project.org')
if (!requireNamespace("rjson", quietly = TRUE))
  install.packages("rjson",repos='http://cran.us.r-project.org')

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(rjson))


scaleFUN <- function(x) sprintf("%.0f", x)

jsonfile<-fromJSON(file = args[1])
kmerlen<-as.character(args[3])
jvals<-as.numeric(paste(unlist(jsonfile)))
jnames<-as.numeric(names(jsonfile))

df<-data.frame(x=jnames,y=jvals,stringsAsFactors = FALSE)

p<-ggplot(df, aes(x,y,group=1)) + geom_point(col="darkred") + geom_line(col="darkred")+ theme_bw()+ labs(x=paste0(kmerlen,'-mer frequency'), y=paste0('# distinct ', kmerlen,'-mers')) + ggtitle(paste0(kmerlen,'-mers spectra')) + theme(plot.title=element_text(hjust=0.5)) + scale_x_continuous(labels=scaleFUN)
ggsave(p, filename=args[2], device="pdf", height=10, width=15)