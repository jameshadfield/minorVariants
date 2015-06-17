#############################################################################################################
# This script takes output from minority_resistance relating sequences to published alleles                 #
# and the ratio of relevent alleles at these positions.                                                     #
#############################################################################################################

rm(list=ls())
# LIBRARIES:
library(ggplot2)
library(reshape2)
# library(RColorBrewer)
library(hash)
library(argparser)

options(error=traceback)
p <- arg.parser("Plot allelic variation over the genome");
p <- add.argument(p, "--seq", help="name of the sequence tab file (python output)")
p <- add.argument(p, "--ref", help="reference tab file (contains genes)")
p <- add.argument(p, "--wd", help="working directory", default=".")
p <- add.argument(p, "--prefix", help="prefix for output")
p <- add.argument(p, "--levels", help="[FLAG] display potential levels", flag=TRUE)
p <- add.argument(p, "--meta", help="additional metadata file")

if (interactive()) {
  argv <- parse.args(p, c("--seq", "./allCTGDmpileupsGZ/K_UK583237.genes.tab.gz","--ref","./het_ref/mutant.v3.tab","--wd","~/Desktop/","--levels","--meta","./het_ref/additional_metadata.tab"));
#   argv <- parse.args(p, c("--seq", "./allCTGDmpileups/L1_L942.genes.tab","--ref","./het_ref/mutant.v3.tab","--wd","~/Desktop/","--levels","--meta","./het_ref/additional_metadata.tab"));
} else {
  argv <- parse.args(p, argv = commandArgs(trailingOnly = TRUE))
}
if ( ! exists("argv") || is.na(argv$seq) || is.na(argv$seq) || is.na(argv$ref)) print(p) & stop()

setwd(argv$wd)
meta <- hash()

#### parse alleles tab file and use this information to shift the "co=ords in gene" to "co-ords in genome" via the offset function\
alleles <- read.table(argv$ref,header=T,sep="",comment.char='', stringsAsFactors=F)
alleles <- alleles[ ! grepl("^#",alleles[,1],), ] 
alleles$offset <- 0
for (idx in seq(1:nrow(alleles))) { alleles[idx,"offset"] <- min(as.numeric(unlist(regmatches(alleles[idx,"location"], gregexpr("[0-9]+", alleles[idx,"location"], perl=TRUE))))) }

#### parse the main data and add in "gpos" which is genomic co-ordinates
if (grepl(".gz$",argv$seq)) {
  raw <- read.table(gzfile(argv$seq),header=T,sep="\t",comment.char='', stringsAsFactors=F)
} else raw <- read.table(argv$seq,header=T,sep="\t",comment.char='', stringsAsFactors=F)
raw$sequence <- factor(raw$sequence)
raw$position <- as.integer(raw$position)
raw$gpos <- 0
h <- hash(alleles[,"geneName"],alleles[,"offset"])
raw$gpos <- apply(raw,1,function(x) h[[ x[[2]] ]] )
raw$gpos <- raw$gpos+raw$position
#### unique(raw$mutation) --> "depth" "WT"    "fixed" "other"

### ensure we are only dealing with one sample here (cannot handle >1)
if ( length(levels(raw$sequence))!=1 ) stop("Can only deal with one sample at a time! FATAL")
meta$seqname <- levels(raw$sequence)[1]

#### split the raw data into two sections: raw.sas and raw.cov
raw.sas <- subset(raw,mutation!="depth")
raw.cov <- subset(raw,mutation=="depth")
raw.sas$type <- "sas" ## for facets
raw.cov$type <- "cov"
meta$mean.depth <- mean(raw.cov$frac)
meta$SNPcallThresh <- 0.95
meta$numSNPs <- nrow(subset(raw,mutation=="fixed" & frac>=meta$SNPcallThresh)) ##### dubious

#### for the sas set, we want to convert to minor -- i.e. map [0,1]->[0,0.5] -- and remove things below some epsilon level (e.g. 1%)
#### this is as we don't know anything about the relatedness of the reference. What follows is reference independent!
epsilon = 0.01
makeAlleleMinor <- function(x) {
  #   x == sequence   name position mutation   frac gpos
  if (as.numeric(x[[5]])<=1 && as.numeric(x[[5]])>0.5) 1 - as.numeric(x[[5]])
  else as.numeric(x[[5]])
}
raw.sas$frac <- apply(raw.sas,1,makeAlleleMinor)
raw.sas <- subset(raw.sas,frac>=epsilon)
if (nrow(raw.sas)==0) stop("No plotting as there are *no* segregating alleles in the whole genome!")

##### for the depth (coverage) set we want a sliding window approach!
slide=500
window=1000
calcWindowAverage <- function(x) {
  #   cat(x[[1]],"\t")
  # x is a vector, specifically a row from the df apply is called from
  as.integer(mean(subset(raw.cov,gpos>=(x[[1]]-(window/2)) & gpos<=(x[[1]]+(window/2)))$frac)) # sloooooow
}
cov.df <- data.frame(x=seq(slide,max(raw.cov$gpos)-slide,by=slide),y=0)
cov.df$y <- apply(cov.df,1,calcWindowAverage)
cov.df$type <- "cov"

#### additional metadata
meta$additional <- ""
if ( ! is.na(argv$meta) ) {
  tmp <- read.table(argv$meta,header=F,sep="\t",comment.char='', stringsAsFactors=F)
  if ( length(levels(raw$sequence))==1 ) {
    meta$additional <- tmp[tmp[,1]==levels(raw$sequence)[[1]],2]
  }
  rm(tmp)
}

####### in raw.sas we have a number of "peaks" in the manhatten plot. Identify which genes these are in
####### and save this information into the regions data.frame
spike_threshold <- 0.25
regions <- data.frame(rName=character(),x1=integer(),x2=integer(),stringsAsFactors=FALSE)
IterBinSearch <- function(A, value) {
  low = 0
  high = nrow(A)
  while ( low <= high ) {
    mid <- floor((low + high)/2)
    coords <- as.numeric(unlist(regmatches(A[mid,"location"], gregexpr("[0-9]+", A[mid,"location"], perl=TRUE))))
    #     cat("value",value,"mid",mid,"coords",coords,"\n")
    if (coords[1] <= value && coords[2] >= value)
      return (mid)
    else if ( coords[1] > value ) #look in the lower half
      high <- mid - 1
    else if ( coords[2] < value ) #look in the upper half
      low <- mid + 1
    else
      break
  }
  #   cat("argh","low",low,"high",high,"mid",mid)
  #   return(A[mid,"geneName"])
  return (NULL)
}
for (idx in seq(1,nrow(raw.sas))) { # each sas variant
  if (raw.sas[idx,"frac"]>spike_threshold) {
    j <- IterBinSearch(alleles,raw.sas[idx,"gpos"]) # what gene this point is "in" using binary search
    if (! is.null(j)) {
#       cat("peak in ",alleles[j,"geneName"],"\n")
      if (! alleles[j,"geneName"] %in% regions[,"rName"]) { # has this gene already been "found"? 
        coords <- as.integer(unlist(regmatches(alleles[j,"location"], gregexpr("[0-9]+", alleles[j,"location"], perl=TRUE))))
        regions[nrow(regions)+1,]  <- c(alleles[j,"geneName"],coords[1],coords[2])
      }
    }
  }
}
regions$x1 <- as.integer(regions$x1)
regions$x2 <- as.integer(regions$x2)
regions$type <- "sas" #so that it goes on the right facet

#### we want to annotate some stuff on the graphs -> data.frame needed (can't annotate per facet)
meta$annstring <- paste("Sequence: ",meta$seqname,"\nMean coverage: ",as.integer(meta$mean.depth),"\nNum called SNPs (>",as.integer(meta$SNPcallThresh*100),"%): ",meta$numSNPs,"\nNum segregating sites (>",as.integer(epsilon*100),"%): ",nrow(raw.sas),"\n","sample details: ",meta$additional,sep="")
meta.df <- data.frame(y=0.49,x=10000,txt=meta$annstring,type="sas")

### some ggplot color stuff. http://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
raw.sas$mutation  <- factor(raw.sas$mutation, levels = c("fixed","WT","other")) # specify order as colours come out red, gree, blue (if n=3)


#### we are now ready to call ggplot
GG <- ggplot() + facet_wrap(~type,ncol=1,scales="free_y")
GG <- GG + geom_point(data=raw.sas,size=2,aes(x=gpos,y=frac,colour=mutation))
GG <- GG + coord_cartesian(xlim=c(0,max(c(alleles$offset,raw.cov$gpos,raw.sas$gpos)))) ## no clean way to do y lims with facets :(  ylim=c(-0.1,0.5))
GG <- GG + geom_text(data=meta.df,aes(x=x,y=y,label=txt,hjust=0,vjust=1))
GG <- GG + geom_text(data=regions,aes(x=x1+(x2-x1)/2,y=-0.05,vjust=1,label=rName,angle=45))
GG <- GG + geom_rect(data=regions, aes(x=x1, y=0, xmin = x1, xmax = x2, ymin = -0.03, ymax = 0, vjust=0), fill = "black")
GG <- GG + geom_line(data=cov.df,aes(x=x,y=y))

### whether we save or display the plot depends on whether it is being called as a script or not!
if (interactive()) {
  GG
} else {
  ggsave(filename=paste(argv$prefix,".SAS.genome.pdf",sep=""), plot = GG, scale = 1, width = 18, height = 12, units = "in", dpi = 300)
  cat("allele plotting along genome finished successfully\n")
}



##################################################################################################################


# 
# ###### ordered on allele frequency --> can we see any levels?
# if (argv$levels && nrow(raw)>50) {
#   newdata <- raw[order(raw$frac),]
#   newdata$neworder <- 0
#   h <- hash(levels(newdata$sequence),rep(0,length(levels(newdata$sequence))))
#   for (idx in seq(1,nrow(newdata))) {
#     seqName <- as.character(newdata[idx,"sequence"])
#     h[[ seqName ]] <- h[[ seqName ]] + 1
#     newdata[idx,"neworder"] <- h[[ seqName ]]                                    
#   }
#   GG2 <- ggplot(data=newdata,aes(x=neworder, y=frac,colour=mutation)) + geom_point() + facet_wrap(~sequence, scales="fixed", ncol=1)
#   
#   if (interactive()) {
#     GG2
#   } else {
#     ggsave(filename=paste(argv$prefix,".SAS.levels.pdf",sep=""), plot = GG2, scale = 1, width = 15, height = 8, units = "in", dpi = 300)
#     cat("alleles ordered by frequency finished successfully\n")
#   }
#   
# }
# 
# 
# # histogram of allele frequencies
# ggplot(data=subset(raw,mutation!="depth" & mutation!="WT" & frac>epsilon),aes(x=frac)) + geom_histogram(stat="bin", binwidth=0.01)
# 
# 
# # histogram of allele frequencies mapped to [0,0.5]
# ggplot(data=raw.sas,aes(x=frac)) + geom_histogram(stat="bin", binwidth=0.005)
# 
# 
# 
# 
# 
# 









