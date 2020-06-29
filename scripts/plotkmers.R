#!/usr/bin/env Rscript

if (!requireNamespace("optparse", quietly = TRUE))
  install.packages("optparse",repos='http://cran.us.r-project.org')
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2",repos='http://cran.us.r-project.org')
if (!requireNamespace("rjson", quietly = TRUE))
  install.packages("rjson",repos='http://cran.us.r-project.org')

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(rjson))

scaleFUN <- function(x) sprintf("%.0f", x)

option_list = list(
  make_option(c("-j", "--json"), action="store", type="character", help=".one or more (comma-separated) .json.gz file from paths kmers [required]"),
  make_option(c("-k", "--kmers"), action="store", default=61, type="numeric", help="k-mers size [61]")
  )

opt = parse_args(OptionParser(option_list=option_list))

kmerlen<-opt$kmers
paths<-unlist(strsplit(opt$json, ','))
pathslist<-list()

for (i in 1:length(paths)) {
	
	p<-file.path(paths[i])
	jsonfile<-fromJSON(file = p)
	jvals<-as.numeric(paste(unlist(jsonfile)))
	jnames<-as.numeric(names(jsonfile))
	name<-basename(p)
	pathslist[[i]]<-cbind(jvals,jnames,name)
	s<-sum(jvals*jnames)
	message("Found ",s, " kmers in ", p)
	peak<-max(jvals)
	indpeak<-which.max(jvals)
	maxkmer<-jnames[indpeak]
	message("Peak is ", maxkmer, ":", peak)
	estimate<-s/maxkmer
	message("Estimated transcriptome size is ", estimate)

}

df<-data.frame(do.call(rbind,pathslist),stringsAsFactors = FALSE)

p<-ggplot(df, aes(x=as.numeric(jvals),y=as.numeric(jnames),group=name,color=name)) + geom_point() + geom_line()+ theme_bw()+ labs(x=paste0(kmerlen,'-mer frequency'), y=paste0('# distinct ', kmerlen,'-mers')) + ggtitle(paste0(kmerlen,'-mers spectra')) + theme(plot.title=element_text(hjust=0.5), legend.title=element_blank()) + scale_color_brewer(palette="Dark2") + scale_x_continuous(labels=scaleFUN)
ggsave(p, filename="kmers.pdf", device="pdf", height=10, width=15)
