# Essentials ----
# Load required packages 
library(tidyverse)
library(edgeR)

# Set working directory 
setwd("D:/AFO/Data/")



# Count data ----
#read count data
counts <- read.delim("GLUCOLD/GLUCOLD - mRNA Gene Expression (RNA Seq)/GLUCOLDhumangeneexpression.txt", 
                     sep = ";", header = T, row.names = 1)

#Filter samples available at baseline
counts_filtered <- counts %>% 
  select(ends_with("_1")) %>% 
  rename_with(~gsub("_1|X", "", .x))



# Patient data ----
#load patient data
patient_data <- haven::read_sav("GLUCOLD/GLUCOLD - Patients/2015-11-24 GLUCOLD database (totaal, 114pt).sav")

#Select required columns
patient_data <- patient_data %>% 
  select(idnr, treatment, geslacht, lftv1, rookjaar, rooknuv1, packyearv1, 
         ccqtotalv1, ccqtotalv3, rvtlcpredv1, rvtlcpredv3, fevnapredv1, fevnapredv3) 

#Filter for patiens that have RNA-seq available at baseline, no NA in any column and have received no placebo
patient_data <- patient_data %>% 
  filter(idnr %in% as.numeric(colnames(counts_filtered))) %>% 
  filter(treatment != 1) %>% 
  drop_na()

#Add a column with the change in lung function after 6 months
patient_data <- patient_data %>% 
  mutate(d_fev_pred = fevnapredv3 - fevnapredv1)

#Final selection of patients 
counts_filtered <- counts_filtered %>% 
  select(which(colnames(counts_filtered) %in% patient_data$idnr))




# RNA-seq data analysis ----
#Design matrix
d_fev_pred <- patient_data$d_fev_pred

mm <- model.matrix(~d_fev_pred)
# mm <- model.matrix(~factor(patient_data$geslacht) + patient_data$lftv1 + d_fev_pred)

row.names(mm) <- patient_data$idnr

#Create a DGEList object
d0 <- DGEList(counts_filtered)

#Normalise counts
d1 <- calcNormFactors(d0)

#estimate dispsion
d2 <- estimateDisp(d1, mm, robust=TRUE)

#
fit <- glmFit(d2, mm)
lrt <- glmLRT(fit)

sign_genes <- topTags(lrt, n=10)
topTags(lrt, n=10)
sum(sign_genes$table$FDR<0.05)


  

