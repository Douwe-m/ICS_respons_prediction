# Load required packages 
library(tidyverse)
library(limma)

#Load Microarray data
expression_data <- read.delim(
  "D:/AFO/Data/GLUCOLD/GLUCOLD - mRNA Gene Expression (Array)/log2glucold221rma2.txt",
  sep = "\t", header = T)

#Gene names
gene_names <- read.delim(
  "D:/AFO/Data/GLUCOLD/GLUCOLD - mRNA Gene Expression (Array)/genenames.txt", 
  sep = "\t", header = T)

#patient data
patient_data <- read.delim(
  "D:/AFO/Data/GLUCOLD/GLUCOLD - mRNA Gene Expression (Array)/s221version3.txt", 
  sep = "\t", header = T) 