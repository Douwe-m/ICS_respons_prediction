# Load required packages 
library(tidyverse)
library(edgeR)
library(GSVA)
library(psych)
library(ggpubr)
library(ggplotify)
library(ggrepel)
library(ggtext)
library(RColorBrewer)
library(pheatmap)
library(patchwork)

data_dir <- "C:/Users/douwe/Documents/Data"
# data_dir <- "D:/AFO/Data"


#read count data
counts <- read.delim(
  file.path(data_dir, "GLUCOLD/GLUCOLD - mRNA Gene Expression (RNA Seq)/GLUCOLDhumangeneexpression.txt"), 
  sep = ";", header = T, row.names = 1)

#Gene symbols
gene_symbols <- read.delim(
  file.path(data_dir, "GLUCOLD/GLUCOLD - mRNA Gene Expression (RNA Seq)/genes_info_hg37.csv"), 
  sep = ",", header = T)

#load patient data
patient_data_raw <- haven::read_sav(
  file.path(data_dir, "GLUCOLD/GLUCOLD - Patients/2015-11-24 GLUCOLD database (totaal, 114pt).sav"))

# #Load normalized counts
# if (file.exists("RNA-seq analysis/output/normalised_counts.csv")){
#   norm_counts1 <- read.table("RNA-seq analysis/output/normalised_counts.csv",
#                             sep = ",",
#                             header = T)
# }

#Load results from DGE analysis
if (file.exists("RNA-seq analysis/output/DEG_d_ccq_v3_baseline.csv")){
  DEG_d_ccq <- read.table("RNA-seq analysis/output/DEG_d_ccq_v3_baseline.csv", 
                          sep = ",", 
                          header = T)
}

if (file.exists("RNA-seq analysis/output/DEG_d_fev_pred_v3_baseline.csv")){
  DEG_d_fev <- read.table("RNA-seq analysis/output/DEG_d_fev_pred_v3_baseline.csv", 
                          sep = ",", 
                          header = T)
}

if (file.exists("RNA-seq analysis/output/DEG_d_rvtlc_v3_baseline.csv")){
  DEG_d_rvtlc <- read.table("RNA-seq analysis/output/DEG_d_rvtlc_v3_baseline.csv", 
                          sep = ",", 
                          header = T)
}
