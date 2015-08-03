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
p <- add.argument(p, "--genome", help="[FLAG] data coming in is from the genome (i.e. gene information is absent)", flag=TRUE)
p <- add.argument(p, "--meta", help="additional metadata file")

if (interactive()) {
  argv <- parse.args(p, c("--seq", "./F_UK583468.genes.tab.gz","--ref","./het_ref/mutant.v3.tab","--wd","~/Desktop/","--genome","--meta","./het_ref/additional_metadata.tab"));
} else {
  argv <- parse.args(p, argv = commandArgs(trailingOnly = TRUE))
}
if ( ! exists("argv") || is.na(argv$seq) || is.na(argv$seq) || is.na(argv$ref)) print(p) & stop()

setwd(argv$wd)
meta <- hash()

#### parse alleles tab file and get minimum and maximum co-ords for each gene
alleles <- read.table(argv$ref,header=T,sep="",comment.char='', stringsAsFactors=F)
alleles <- alleles[ ! grepl("^#",alleles[,1],), ]
alleles$x1 <- 0 ## minimum value of a gene region (e.g. smallest number)
alleles$x2 <- 0
for (idx in seq(1:nrow(alleles))) {
  tmp <- as.numeric(unlist(regmatches(alleles[idx,"location"], gregexpr("[0-9]+", alleles[idx,"location"], perl=TRUE))))
  alleles[idx,"x1"] <- min(tmp)
  alleles[idx,"x2"] <- max(tmp)
}


#### parse the main data
if (grepl(".gz$",argv$seq)) {
  raw <- read.table(gzfile(argv$seq),header=T,sep="\t",comment.char='', stringsAsFactors=F)
} else raw <- read.table(argv$seq,header=T,sep="\t",comment.char='', stringsAsFactors=F)
colnames(raw) <- c("sequence", "genename", "position", "WT", 'ALT', 'AT>TA', 'AT>CG', 'GC>CG', 'GC>TA', 'AT>GC', 'GC>AT', 'depth')
raw$sequence <- factor(raw$sequence)
raw$position <- as.integer(raw$position)

#### create the variable "gpos" which is the genome position. This depends on the incoming data format!
if (argv$genome) {
  raw$gpos <- raw$position
} else {
  raw$gpos <- 0
  h <- hash(alleles[,"geneName"],alleles[,"x1"])
  raw$gpos <- apply(raw,1,function(x) h[[ x[[2]] ]] )
  raw$gpos <- raw$gpos+raw$position
}

### ensure we are only dealing with one sample here (cannot handle >1)
if ( length(levels(raw$sequence))!=1 ) stop("Can only deal with one sample at a time! FATAL")
meta$seqname <- levels(raw$sequence)[1]
meta$mean.depth <- as.integer(mean(raw$depth))
meta$SNPcallThresh <- 0.95
meta$numSNPs <- sum(raw$ALT > meta$SNPcallThresh) ### different to T81 scripts.

#### additional metadata
meta$additional <- ""
if ( ! is.na(argv$meta) ) {
  tmp <- read.table(argv$meta,header=F,sep="\t",comment.char='', stringsAsFactors=F)
  if ( length(levels(raw$sequence))==1 ) {
    meta$additional <- tmp[tmp[,1]==levels(raw$sequence)[[1]],2]
  }
  rm(tmp)
}

##### create a data.frame to have the depth information via sliding window
slide=5000
window=10000
calcWindowAverage <- function(x) {
  #   cat(x[[1]],"\t")
  # x is a vector, specifically a row from the df apply is called from
  as.integer(mean(raw[raw$gpos>=(x[[1]]-(window/2)) & raw$gpos<=(x[[1]]+(window/2)),"depth"]))
}
cov.df <- data.frame(x=seq(slide,max(raw$gpos)-slide,by=slide),y=0)
cov.df$y <- apply(cov.df,1,calcWindowAverage)
cov.df$type <- "cov"
cov.df <- cov.df[! is.na(cov.df$y),]# remove missing values

###### say we want to look at transitions and transversions....
main <- raw[,c("genename","gpos")]
main$transitions <- apply(raw[,c('AT>GC', 'GC>AT')],1,sum)
main$transversions <- apply(raw[,c('AT>TA', 'AT>CG', 'GC>CG', 'GC>TA')],1,sum)


##### convert to minor -- i.e. map [0,1]->[0,0.5] -- and remove things below some epsilon level for plotting purposes (e.g. 1%)
epsilon = 0.01
makeAlleleMinor <- function(x,colnum) {
  #   x == sequence   name position mutation   frac gpos
  if (as.numeric(x[[colnum]])<=1 && as.numeric(x[[colnum]])>0.5) 1 - as.numeric(x[[colnum]])
  else as.numeric(x[[colnum]])
}
main[,3] <- apply(main,1,makeAlleleMinor,colnum=3)
main[,4] <- apply(main,1,makeAlleleMinor,colnum=4)
main <- melt(main,id.vars=c("genename","gpos"))
main <- subset(main, value>=epsilon )
if (nrow(main)==0) stop("No plotting as there are *no* segregating alleles in the whole genome!")
meta$numSAS <- length(unique(main$gpos))
main$type = "sas"
meta$changes <- paste( "K: ", sprintf("%.2f", sum(main$variable=="transitions") / sum(main$variable=="transversions") ), sep="")

#### if argv$genome then we don't have the gene information -- just the bases in the genome. So add in *the ability to get this information* via a binary search.
if (argv$genome) {
  identify_gene_from_genome_position_via_binary_search <- function(A,value) {
    high_idx = nrow(A)
    low_idx = 0
    while (low_idx<=high_idx) {
      mid_idx  <- ceiling((low_idx+high_idx)/2)
      if (value < A[mid_idx,"x1"]) { # enter lower halp
        high_idx <- mid_idx - 1
      } else if  (value > A[mid_idx,"x2"]) { # enter upper half
        low_idx <- mid_idx + 1
      } else if (value >= A[mid_idx,"x1"] && value <= A[mid_idx,"x2"]) { # hit
        return(A[mid_idx,"geneName"])
      } else {
#         cat("1intergenic gene? mid position is",mid_idx,"value is",value,"df row is",as.character(A[mid_idx,]),"\n")
        return("intergenic")
      }
    }
#     cat("2intergenic gene? mid position is",mid_idx,"value is",value,"df row is",as.character(A[mid_idx,]),"\n")
    return("intergenic")
  }
}


##### from the manhatten-like plot, we want all the "spikes" a.k.a. collections of points over some value to have their corresponding gene displayed...
spike_threshold <- 0.25
regions <- data.frame(rName=character(),x1=integer(),x2=integer(),stringsAsFactors=FALSE)
for (idx in seq(1,nrow(main))) { # each SaS variant
  if (main[idx,"value"]>spike_threshold) {
    if (argv$genome) {
      geneName <- identify_gene_from_genome_position_via_binary_search(alleles,main[idx,"gpos"])
    } else {
      geneName <- main[idx,"genename"]
    }
    if (! geneName %in% regions[,"rName"]) { # has this gene already been "found"?
      if (geneName=="intergenic") {
        regions[nrow(regions)+1,]  <- c(geneName,main[idx,"gpos"],main[idx,"gpos"])
      } else {
        regions[nrow(regions)+1,]  <- c(geneName,alleles[alleles$geneName==geneName,"x1"],alleles[alleles$geneName==geneName,"x2"])
      }
    }
  }
}
regions$x1 <- as.integer(regions$x1)
regions$x2 <- as.integer(regions$x2)
regions$type <- "sas" #so that it goes on the right facet


#### we want to annotate some stuff on the graphs -> data.frame needed (can't annotate per facet)
meta$annstring <- paste("Sequence: ",meta$seqname,"\nMean coverage: ",as.integer(meta$mean.depth),"\nNum called SNPs (>",as.integer(meta$SNPcallThresh*100),"%): ",meta$numSNPs,"\nNum segregating sites (>",as.integer(epsilon*100),"%): ",meta$numSAS,"\n","sample details: ",meta$additional,"\n",meta$changes,sep="")
meta.df <- data.frame(y=0.49,x=10000,txt=meta$annstring,type="sas")

### some ggplot color stuff. http://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length=n+1)
#   hcl(h=hues, l=65, c=100)[1:n]
# }
# raw.sas$mutation  <- factor(raw.sas$mutation, levels = c("fixed","WT","other")) # specify order as colours come out red, gree, blue (if n=3)
#


#### we are now ready to call ggplot
GG <- ggplot() + facet_wrap(~type,ncol=1,scales="free_y")
GG <- GG + geom_point(data=main,size=2,aes(x=gpos,y=value,colour=variable))
# GG <- GG + coord_cartesian(xlim=c(0,max(c(alleles$offset,raw.cov$gpos,raw.sas$gpos)))) ## no clean way to do y lims with facets :(  ylim=c(-0.1,0.5))
GG <- GG + geom_text(data=meta.df,aes(x=x,y=y,label=txt,hjust=0,vjust=1))
GG <- GG + geom_text(data=regions,aes(x=x1+(x2-x1)/2,y=-0.05,vjust=1,label=rName,angle=45))
GG <- GG + geom_rect(data=regions, aes(x=x1, y=0, xmin = x1, xmax = x2, ymin = -0.03, ymax = 0, vjust=0), fill = "black")
GG <- GG + geom_line(data=cov.df,aes(x=x,y=y))
GG <- GG + scale_colour_manual(values = c("#41b6c4", "#2c7fb8"))

### whether we save or display the plot depends on whether it is being called as a script or not!
if (interactive()) {
  GG
} else {
  ggsave(filename=paste(argv$prefix,".SAS.genome.pdf",sep=""), plot = GG, scale = 1, width = 18, height = 12, units = "in", dpi = 300)
  cat("allele plotting along genome finished successfully\n")
}


###### allele frequency plot:
epsilon <- 0.01
base.changes.vec <- c('AT>TA', 'AT>CG', 'GC>CG', 'GC>TA', 'AT>GC', 'GC>AT')
AF <- data.frame(x=base.changes.vec,y=0,class=factor(c( rep("transversion",4), rep("transition",2) ), levels=c("transversion","transition")),stringsAsFactors=F)
for (i in AF$x) AF[which(AF$x==i),"y"] <- length(raw[raw[i]>epsilon & raw[i]<1-epsilon,i])
AF$x <- factor(AF$x,levels=base.changes.vec)
total_sites <- sum(AF[,"y"])
for (i in seq(1,nrow(AF))) AF[i,"y"] <- AF[i,"y"]/total_sites
ggplot(data=AF,aes(x=x,y=y))+geom_bar(stat="identity",aes(fill=class))

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









