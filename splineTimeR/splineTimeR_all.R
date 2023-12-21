require(splineTimeR)
require(dplyr)
require(ggplot2)

#------------------------------------------------------------------------------#
setwd("/home/zhiyong/Desktop/ZZZZZZZZZZZ")   #this is the directery this package will save its plots in

#------------------------------------------------------------------------------#
data(TCsimData)   #load example data

pData(TCsimData)   #check the experimental design

#------------------------------------------------------------------------------#
#detecting differential expression genes along time-curse between two condition
diffExprs <- splineDiffExprs(eSetObject=TCsimData, df=3, cutoff.adj.pVal=0.01, reference="T1", intercept=F)   #df=

class(diffExprs)
dim(diffExprs)
head(diffExprs)

#------------------------------------------------------------------------------#
#visualization of differential expression genes along time-curse between two condition
splinePlot(eSetObject = TCsimData, df = 3, reference = "T1", toPlot = c("EEF1G"))