#############################################################################################################
# This script takes output from minority_resistance relating sequences to published alleles                 #
# and the ratio of relevent alleles at these positions.                                                     #
# It can be run interactively (by changing values in the "parameters" section)                              #
# or via Rscript with the following options:                                                                #
#   RScript plot_minor_variants.R working_directory tabfileoutput tabfileinput gene output plot_params      #
#                                        1               2            3         4       5       6           #
#       gene:         gene name (see plot_params) (this *must* be set, even if unused)                      #
#       tabfileout:   output from minority_resistance                                                       #
#       tabfilein:    input to minority_resistance                                                          #
#       output:       file to save (with pdf extension)                                                     #
#       plot_params:  gene(1|0),numColumns,maxYvalue                                                        #
#############################################################################################################

###### P A R A M E T E R S #####
if (length(commandArgs(trailingOnly=TRUE))) {
  interactive  <-  FALSE
  workDir      <-  commandArgs(TRUE)[[1]]
  fname        <-  commandArgs(TRUE)[[2]]
  tabName      <-  commandArgs(TRUE)[[3]]
  geneName     <-  commandArgs(TRUE)[[4]] 
  saveName     <-  commandArgs(TRUE)[[5]]
  plot_params  <-  as.vector(sapply(strsplit(commandArgs(TRUE)[[6]], ",")[[1]],as.numeric))
  geneAnalysis <-  as.logical(plot_params[1])  
  ncol         <-  plot_params[2]
  ylimUpper    <-  plot_params[3]
} else { ## interactive // testing --> this section changes frequently
  interactive  <-  TRUE
  workDir      <-  "/Volumes/user_homes_1b/nfs_j/jh22/tmp/23S_mutations"
  fname        <-  "H.alleles.tab"
  tabName      <-  "23Smutations.tab"
  geneName     <-  ""
  geneAnalysis <-  FALSE
  ncol         <-  2
  ylimUpper    <-  0.2
}

# LIBRARIES:
library(ggplot2)
library(reshape2)
library(RColorBrewer)

########## FUNCTIONS ##########
draw <- function(resistancedf) {
  #### first we tidy up the xlabels
  xLabels <- sort(unique(as.character(resistancedf$position)))
  if (! geneAnalysis) { ## attach drug information to X labels
    ## quickly parse the tabfile to get the bindings of gene -> drug (if applicable)
    drugMap <- read.table(tabName,header=T,sep="",comment.char='', stringsAsFactors=F)[,c(1,10)]
    stack <- vector()
    for (label in xLabels) {
      drug <- Filter(function(x) x!="-",drugMap[drugMap[,1]==strsplit(label, "_")[[1]][[1]],2])[1]
      if (is.na(drug)) {drug = ""}
      stack <- c(stack, drug) 
    }
    xLabels <- paste(xLabels,stack,sep = "\n")
  }
    
  ## ordering for ggplot:
  resistancedf$mutation <- factor(resistancedf$mutation,levels=c("published","other","WT","fixed"), ordered=TRUE)
  resistancedf$position <- factor(resistancedf$position, levels=sort(unique(as.character(resistancedf$position))), ordered=TRUE)
  
  ## ggplot
  GG <-      ggplot(resistancedf,aes(x=position, y=frac, fill = mutation, order=mutation))
  GG <- GG + geom_bar(stat="identity")
  GG <- GG + facet_grid(sequence ~ . )
  # GG <- GG + facet_wrap(~ sequence , ncol=2)
  GG <- GG + theme(axis.text.x=element_text(angle = -45, hjust = 0))
  GG <- GG + ylab("fraction of reads")
  # GG <- GG + scale_fill_manual(values=c("red", "blue", "white"))
  GG <- GG + coord_cartesian(ylim = c(0,ylimUpper)) 
  # GG <- GG + guides(fill = guide_legend(reverse = TRUE))
  GG <- GG + scale_fill_manual(values = rev(brewer.pal(3,"PuRd")), name="AA Mutation")
  # GG <- GG + scale_fill_brewer(palette="YlGnBu")
  GG <- GG + theme(legend.position="bottom")
  if (geneAnalysis) {
    GG <- GG + scale_x_discrete(breaks=NULL)
    GG <- GG + xlab(geneName)
  } else {
    GG <- GG + scale_x_discrete(labels = xLabels) ## scale_x_discrete(breaks=as.integer(seq(from=1,to=length(xLabels),length.out=10)) ),labels=xLabels[as.integer(seq(from=1,to=length(xLabels),length.out=10))])
    GG <- GG + xlab("mutation")
  }
  GG <- GG +  scale_y_continuous(breaks = seq(0,ylimUpper-0.1,0.1))
  
  return(GG)
}

# Multiple plot function
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),ncol = cols, nrow = ceiling(numPlots/cols))  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,layout.pos.col = matchidx$col))
    }
  }
}


setwd(workDir)
resistancedf <- read.table(fname,header=T,sep="\t",comment.char='', stringsAsFactors=F)
resistancedf <- resistancedf[resistancedf$mutation!="WT",] #get rid of wt fraction (else all bars range [0,1])
if (geneAnalysis) {  # restrict to chosen gene name
  resistancedf <- resistancedf[grepl(geneName,resistancedf$position),]
}

### the plot is quite easy to draw by calling draw(someDF) -> GGplot object (e.g. GG <- draw(resistancedf))
### but we need to tweak formatting (e.g. num cols -> multiplot) and if we should save it!
if ( ! interactive ){
  height = 0.66 * length(unique(resistancedf[ , 1])) ## .66 inch per track
  width = 20 #0.33 * length(levels(resistancedf$position)) * ncol
  # fit to A4 page
  if (width > 20) {width = 20} 
  if (height > 11) {height = 11}
  # minimum sizes:
  if (width < 6) {width = 6}
  if (height < 6) {height = 6}
  pdf(saveName, width = width, height = height)
}

print(height)
print(width)

## if asked for we can display multiple plots side by side...
seqs  <- unique(resistancedf$sequence)
npercol  <- ceiling(length(seqs)/ncol)
GGlist <- list()
endVal <- 0
for (i in 1:ncol) {
  startVal<- endVal + 1
  endVal  <- startVal + npercol - 1
  print (c(length(seqs),startVal,endVal))
  print(resistancedf[resistancedf$sequence %in% seqs[startVal:endVal],])
  GG <- draw(resistancedf[resistancedf$sequence %in% seqs[startVal:endVal],])
  GGlist[[i]]  <- GG
}
multiplot(plotlist=GGlist, layout=matrix(seq(1,ncol), nrow=1, byrow=TRUE))

if ( ! interactive ){ # this saves the pdf
  dev.off()
}





# ## e.g. how many strains have a "published" mutation >5% in rplD (azithro)?
# length(unique(resistancedf[resistancedf$position=="rplD_66" & resistancedf$mutation=="published" & resistancedf$frac>=0.05 , 1]))
# ## e.g. how many strains have a "other" mutation >5% in rplD (azithro)?
# length(unique(resistancedf[resistancedf$position=="rplD_66" & resistancedf$mutation=="other" & resistancedf$frac>=0.05 , 1]))
# 
# ## e.g. how many strains have a "other" mutation >5% in rpoB ?
# length(unique(resistancedf[grepl('rpoB',resistancedf$position) & resistancedf$mutation=="other" & resistancedf$frac>=0.05 , 1]))
# ## e.g. how many strains have a "published" mutation >5% in rpoB ?
# length(unique(resistancedf[grepl('rpoB',resistancedf$position) & resistancedf$mutation=="published" & resistancedf$frac>=0.05 , 1]))






