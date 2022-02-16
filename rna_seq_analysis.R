# Essentials ----
# Load required packages 
library(tidyverse)

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

# Patient characteristics ----
# Calculate patient characteristics
#Gender
patient_data %>% 
  group_by(geslacht) %>% 
  count(geslacht)

#Age
mean(patient_data$lftv1)
sd(patient_data$lftv1)
 
#smoking status
patient_data %>% 
   group_by(rooknuv1) %>% 
   count(rooknuv1)
 
#Smoking history
mean(patient_data$packyearv1)
sd(patient_data$packyearv1)
 
#FEV1 % pred baseline
patient_data %>% 
  summarise_at(vars(fevnapredv1, fevnapredv3, ccqtotalv1, ccqtotalv3, rvtlcpredv1, rvtlcpredv3), 
               c(mean = "mean", sd = "sd")) %>% 
  t()













# Other stuff ----
#Select samples at baseline
subset_counts <- counts %>% select(as.character(patient_data$idnr))

d0 <- DGEList(subset_counts)
d1 <- calcNormFactors(d0)

mm <- model.matrix(~patient_data$geslacht + patient_data$lftv1 + patient_data$delta)
row.names(mm) <- patient_data$idnr

#????
y <- estimateDisp(d1, mm, robust=TRUE)
fit <- glmFit(y, mm)
lrt <- glmLRT(fit)
topTags(lrt)
