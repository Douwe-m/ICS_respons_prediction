#Load required packages
library(edgeR)

#Set working directory
setwd("D:/AFO/Data/")

counts <- read.table("GLUCOLD/GLUCOLD - mRNA Gene Expression (RNA Seq)/GLUCOLDhumangeneexpression.txt", header = T, sep = ";", row.names = 1)
