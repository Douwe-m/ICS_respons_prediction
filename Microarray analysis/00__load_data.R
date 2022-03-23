# Load required packages 
library(tidyverse)
library(limma)
library(here)
library(GSVA)

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


#Load results from DGE analysis
if (file.exists(here("Microarray analysis", "output", "DEG_d_ccq_v3_baseline.csv"))){
  DEG_d_ccq <- read.table(here("Microarray analysis", "output", "DEG_d_ccq_v3_baseline.csv"),
                          sep = ",",
                          header = T)
}
if (file.exists(here("Microarray analysis", "output", "DEG_d_fev_pred_v3_baseline.csv"))){
  DEG_d_fev <- read.table(here("Microarray analysis", "output", "DEG_d_fev_pred_v3_baseline.csv"),
                          sep = ",",
                          header = T)
}
if (file.exists(here("Microarray analysis", "output", "DEG_d_rvtlc_v3_baseline.csv"))){
  DEG_d_rvtlc <- read.table(here("Microarray analysis", "output", "DEG_d_rvtlc_v3_baseline.csv"),
                          sep = ",",
                          header = T)
}
