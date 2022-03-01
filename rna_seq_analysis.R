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

#Gene symbols
gene_symbols <- read.delim("D:/AFO/Data/GLUCOLD/GLUCOLD - mRNA Gene Expression (RNA Seq)/genes_info_hg37.csv", sep = ",", header = T)

# Patient data ----
#load patient data
patient_data <- haven::read_sav("GLUCOLD/GLUCOLD - Patients/2015-11-24 GLUCOLD database (totaal, 114pt).sav")

#Select required columns
patient_data <- patient_data %>% 
  select(idnr, treatment, geslacht, lftv1, rookjaar, rooknuv1, packyearv1, 
         ccqtotalv1, ccqtotalv3, rvtlcpredv1, rvtlcpredv3, fevnapredv1, fevnapredv3) 

#Filter for patients that have RNA-seq available at baseline, no NA in any column and have received no placebo
patient_data <- patient_data %>% 
  filter(idnr %in% as.numeric(colnames(counts_filtered))) %>% 
  filter(treatment != 1) %>% 
  drop_na()

#Add a column with the change in lung function after 6 months
patient_data <- patient_data %>% 
  mutate(d_fev_pred = fevnapredv3 - fevnapredv1) %>% 
  mutate(d_ccq = ccqtotalv3 - ccqtotalv1) %>% 
  mutate(d_rvtlc = rvtlcpredv3 - rvtlcpredv1)


#Final selection of patients
counts_filtered <- counts_filtered %>%
  select(which(colnames(counts_filtered) %in% patient_data$idnr))

top_genes <- list()

# Analyse RNA-seq data ----
for (i in select(patient_data, starts_with("d_"))) {
  #Design matrix
  mm <- model.matrix(~i)
  row.names(mm) <- patient_data$idnr

  #Create a DGEList object
  d0 <- DGEList(counts_filtered)
  keep <- filterByExpr(d0)
  d0 <- d0[keep, , keep.lib.sizes=FALSE]
  
  #Normalize counts
  d1 <- calcNormFactors(d0)
  
  #estimate dispersion
  d2 <- estimateDisp(d1, mm, robust=TRUE)
  
  fit <- glmQLFit(d2, mm)
  fit <- glmQLFTest(fit)
  
  tt <- topTags(fit, n = nrow(d2), adjust.method = "BH")
  
  top_tt <- tt$table %>%
    rownames_to_column("ensembl_gene_id") %>%
    left_join(y = gene_symbols %>% select(ensembl_gene_id, hgnc_symbol), by = "ensembl_gene_id") %>% 
    distinct()
  
  top_genes <- c(list(top_tt), top_genes)
}

