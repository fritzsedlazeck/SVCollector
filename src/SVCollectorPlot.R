#!/usr/bin/env Rscript
  
COL_GREEDY="red"
COL_TOPN="orange"
COL_RANDOM="black"

COL_BGCOLOR="light grey"
COL_AXIS="white"

args<-commandArgs(TRUE)

if (length(args) < 3) {
  cat("USAGE: SVCollector.R svcollector.greedy svcollector.topN svcollector.random out.png\n")
} else {
  greedyf <- args[[1]]
  topNf   <- args[[2]]
  randomf <- args[[3]]
  outf    <- args[[4]]

  cat("  loading greedy:", greedyf, "\n")
  greedy <- read.table(greedyf, skip=1)

  cat("  loading topN:", topNf, "\n")
  topN   <- read.table(topNf,   skip=1)

  cat("  loading random:", randomf, "\n");

  randomdir    = dirname(randomf)
  randomprefix = basename(randomf)
  randomscan   = paste(randomprefix, "\\.", sep="")
  randomfiles  = list.files(path=randomdir, pattern=randomscan)

  idx = 0

  for (f in randomfiles) {
    fp = paste(randomdir,"/",f, sep="")
    cat("  loading random file: ", fp, "\n")

    d <- read.table(fp, skip=1)

    if (idx == 0)
    {
      random = d[[4]]
      idx = idx + 1
    }
    else
    {
      random = cbind(random, d[[4]])
    }
  }

  numsamples = max(dim(greedy))

  cat("ploting to:", outf, "\n")

  png(outf, height=1000, width=1000)

  plot(greedy[[4]],  typ="n", ylim=c(0,1), 
       ylab="Cummulative Fraction of SVs", xlab="Num Samples", main=paste("SVCollector for", substr(randomprefix, 1, nchar(randomprefix)-7)))

  rect(-10, -10, numsamples*1.5, 2, col=COL_BGCOLOR)
  abline(h=c(0,.2,.4,.6,.8,1), col=COL_AXIS)
  xspace=floor(numsamples/5)
  abline(v=seq(0,numsamples,xspace), col=COL_AXIS)

  boxplot(t(random), add=TRUE)

  lines(greedy[[4]],  typ="l", col=COL_GREEDY, lwd=3)

  lines(topN[[4]],    typ="l", col=COL_TOPN, lwd=3)

  legend("topleft", inset=0.05, 
         legend=c("greedy", "topN", "random"), 
         fill=c(COL_GREEDY, COL_TOPN, COL_RANDOM))

  trapmsg <- dev.off()
}
