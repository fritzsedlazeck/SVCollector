#!/usr/bin/env Rscript

args<-commandArgs(TRUE)

if (length(args) < 2) {
  cat("USAGE: plotfit.R svcollector.greedy out.png\n")
} else {
  dataf <- args[[1]]
  outf  <- args[[2]]

  cat ("loading: ", dataf, "\n")
  t <- read.table(dataf, skip=1)
  y = t[[4]]
  x <- seq(length(y))

  cat(head(x), "\n")
  cat(head(y),"\n")
  data <- data.frame(x,y)


  cat ("plotting to:" , outf, "\n")
  png(outf)

  plot(x,y)
  # starting value
  d = 0.03 
  
  # slope at top
  tg = 0.012

## initial slope
  ks = 0.045

## distance for initial segment
  ts = 4.7

## Equation 4 from: http://www.pnas.org/content/pnas/102/39/13950.full.pdf
  pn = d + tg * (x-1) + ks * exp(-2/ts) * (1-exp(-(x-1)/ts)) / (1-exp(-1/ts))
  lines(x,pn, col="red")
  

  #model <- nls ( y ~ (d + tg * (x-1) + ks * exp(-2/ts) * (1-exp(-x-1)/ts) / (1-exp(-1/ts))), data, start=list(d=.03, tg=.012, ts=4.5, ks=.05))
  #cat(model)


  msg <- dev.off()
}

