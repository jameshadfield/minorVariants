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

if (interactive()) {
  argv <- parse.args(p, c("--seq", "./allCTGDmpileups/D_NL17.genes.tab","--ref","./het_ref/mutant.v3.tab","--wd","~/Desktop/","--levels"));
} else {
  argv <- parse.args(p, argv = commandArgs(trailingOnly = TRUE))
}
if ( ! exists("argv") || is.na(argv$seq) || is.na(argv$seq) || is.na(argv$ref)) print(p) & stop()

setwd(argv$wd)


#### parse alleles tab file and use this information to shift the "co=ords in gene" to "co-ords in genome" via the offset function\
alleles <- read.table(argv$ref,header=T,sep="",comment.char='', stringsAsFactors=F)
alleles <- alleles[ ! grepl("^#",alleles[,1],), ] 
alleles$offset <- 0
for (idx in seq(1:nrow(alleles))) { alleles[idx,"offset"] <- min(as.numeric(unlist(regmatches(alleles[idx,"location"], gregexpr("[0-9]+", alleles[idx,"location"], perl=TRUE))))) }

raw <- read.table(argv$seq,header=T,sep="\t",comment.char='', stringsAsFactors=F)
raw$sequence <- factor(raw$sequence)
raw$position <- as.integer(raw$position)
raw <- subset(subset(raw,mutation!="WT"),mutation!="depth")
raw$gpos <- 0
h <- hash(alleles[,"geneName"],alleles[,"offset"])
raw$gpos <- apply(raw,1,function(x) h[[ x[[2]] ]] )
raw$gpos <- raw$gpos+raw$position
# for (idx in seq(1:nrow(raw))) { cat(idx,"\t"); raw[idx,"gpos"] <- raw[idx,"position"] + min(alleles[alleles$geneName==raw[idx,"name"],"offset"] )} ### slooow
# for (idx in seq(1:nrow(raw))) { cat(idx,"\t"); raw[idx,"gpos"] <- raw[idx,"position"] + h[[ raw[idx,"name"] ]] }

### plotting the alleles frequencies across the chromosome to see any enriched areas
#   (geom_point as geom_bar is too sparse!)
#   the gene names e.t.c. are shown too (maybe intelligently?)
makeAlleleMinor <- function(x) {
  #   x == sequence   name position mutation   frac gpos
  if (x[[5]]>0.5) 1 - as.numeric(x[[5]])
  else as.numeric(x[[5]])
}
raw$frac <- apply(raw,1,makeAlleleMinor)
epsilon <- 0.01  ## 1%
raw <- subset(raw,frac>=epsilon)
if (nrow(raw)==0) stop("No plotting as there are *no* segregating alleles in the whole genome!")
# for (idx in seq(1,nrow(newdata))) { if (newdata[idx,"frac"]>0.5) {newdata[idx,"frac"] <- 1 - newdata[idx,"frac"] }  } ## slow
GG <- ggplot(data=raw,aes(x=gpos, y=frac)) + geom_point(size=2,aes(colour=mutation))
GG <- GG + facet_wrap(~sequence, scales="fixed", ncol=1)
GG <- GG + coord_cartesian(xlim=c(0,max(c(alleles$offset,raw$gpos))))

##### from the manhatten-like plot, we want all the "spikes" a.k.a. collections of points over some value to have their corresponding gene displayed...
spike_threshold <- 0.25
regions <- data.frame(rName=character(),x1=integer(),x2=integer(),sequence=character(),stringsAsFactors=FALSE)
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
for (idx in seq(1,nrow(raw))) {
  if (raw[idx,"frac"]>spike_threshold) {
    ### find out what gene this point is "in" using binary search
    j <- IterBinSearch(alleles,raw[idx,"gpos"])
    if (! is.null(j)) {
      ### has this gene already been "found"? (in this sequence only -- faceting)
      sequence <- as.character(raw[idx,"sequence"])
      if ( ! alleles[j,"geneName"] %in% regions[regions$sequence==sequence,"rName"] ) {
        ### add into a dataframe (regions, colu mns: gName, x1, x2)
        coords <- as.integer(unlist(regmatches(alleles[j,"location"], gregexpr("[0-9]+", alleles[j,"location"], perl=TRUE))))
        #         cat("SPIKE!","at",raw[idx,"gpos"],"in",alleles[j,"geneName"],alleles[j,"location"],"sequence",as.character(raw[idx,"sequence"]),"\n")
        regions[nrow(regions)+1,] <- c(alleles[j,"geneName"],coords[1],coords[2],as.character(raw[idx,"sequence"]))
      }      
    } 
    #     else {
    #       cat("SPIKE!","at",raw[idx,"gpos"],"not in a gene...\n")
    #     }
  }
}
regions$x1 <- as.integer(regions$x1)
regions$x2 <- as.integer(regions$x2)
regions$sequence <- as.factor(regions$sequence)


GG <- GG + geom_text(data=regions,aes(x=x1+(x2-x1)/2,y=-0.1,label=rName,angle=45)) + geom_rect(data=regions, aes(x=x1, y=0, xmin = x1, xmax = x2, ymin = -0.05, ymax = -0.01), fill = "black")
### whether we save or display the plot depends on whether it is being called as a script or not!
if (interactive()) {
  GG
} else {
  ggsave(filename=paste(argv$prefix,".SAS.genome.pdf",sep=""), plot = GG, scale = 1, width = 15, height = 8, units = "in", dpi = 300)
  cat("allele plotting along genome finished successfully\n")
}



###### ordered on allele frequency --> can we see any levels?
if (argv$levels && nrow(raw)>50) {
  newdata <- raw[order(raw$frac),]
  newdata$neworder <- 0
  h <- hash(levels(newdata$sequence),rep(0,length(levels(newdata$sequence))))
  for (idx in seq(1,nrow(newdata))) {
    seqName <- as.character(newdata[idx,"sequence"])
    h[[ seqName ]] <- h[[ seqName ]] + 1
    newdata[idx,"neworder"] <- h[[ seqName ]]                                    
  }
  GG2 <- ggplot(data=newdata,aes(x=neworder, y=frac,colour=mutation)) + geom_point() + facet_wrap(~sequence, scales="fixed", ncol=1)
  
  if (interactive()) {
    GG2
  } else {
    ggsave(filename=paste(argv$prefix,".SAS.levels.pdf",sep=""), plot = GG2, scale = 1, width = 15, height = 8, units = "in", dpi = 300)
    cat("alleles ordered by frequency finished successfully\n")
  }
  
}



















