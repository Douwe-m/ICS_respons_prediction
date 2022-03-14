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

#Create an empty list
top_genes <- list()

# Analyse RNA-seq data ----
#Create a DGEList object
d0 <- DGEList(counts_filtered)
keep <- filterByExpr(d0)
d0 <- d0[keep, , keep.lib.sizes = F]

#Normalize counts
d1 <- calcNormFactors(d0)

#Get TMM normalized counts
norm_counts <- cpm(d1, log = T)

for (i in select(patient_data, starts_with("d_"))) {
  #Design matrix
  mm <- model.matrix(~i)
  row.names(mm) <- patient_data$idnr

  #estimate dispersion
  d2 <- estimateDisp(d1, mm, robust = T)
  
  fit <- glmQLFit(d2, mm)
  fit <- glmQLFTest(fit)
  
  tt <- topTags(fit, n = nrow(d2), adjust.method = "BH")
  
  top_tt <- tt$table %>%
    rownames_to_column("ensembl_gene_id") %>%
    left_join(y = gene_symbols %>% select(ensembl_gene_id, hgnc_symbol), by = "ensembl_gene_id") %>% 
    distinct()
  
  top_genes <- c(list(top_tt), top_genes)

}

# Visualizations ----
for (i in top_genes){
  #Check for significant genes
  if (sum(i$FDR <= 0.05) > 0) {
  
    #Select genes with FDR <= 0.05
    sign_genes <- i %>% 
      filter(FDR <= 0.05)
    
    #create df with counts of sign. genes and clinical parameter
    a <- norm_counts %>% 
      as.data.frame() %>% 
      rownames_to_column("geneID") %>% 
      filter(geneID %in% sign_genes$ensembl_gene_id) %>% 
      pivot_longer(2:37, names_to = "patient", values_to = "expr") %>% 
      mutate(patient = as.numeric(patient)) %>% 
      left_join(y = patient_data %>% select(idnr, d_ccq, ccqtotalv1, ccqtotalv3), by = c("patient" = "idnr"))
  }
}

#scatter plot
ggplot(a, aes(y = d_ccq, x = expr, group = factor(geneID))) + 
  geom_point() +
  labs(x = "Log2 TMM values", y = expression(Delta ~ "CCQ Total Score"))+
  facet_wrap(. ~ geneID, ncol = 2)


#boxplot
b <- a %>% 
  pivot_longer(5:6, names_to = "measurement", values_to = "CCQ") %>% 
  mutate(measurement = as.factor(measurement))

levels(b$measurement) <- c("baseline", "6_months") 

b %>% 
  filter(geneID == "ENSG00000133110") %>% 
  ggplot(aes(x = expr, y = CCQ, group = measurement)) +
  geom_point() +
  facet_wrap(. ~ measurement, ncol = 2)

b %>% 
  filter(geneID == "ENSG00000161905") %>% 
  ggplot(aes(x = expr, y = CCQ, group = measurement)) +
  geom_point() +
  facet_wrap(. ~ measurement, ncol = 2)

b %>% 
  filter(geneID == "ENSG00000264940") %>% 
  ggplot(aes(x = expr, y = CCQ, group = measurement)) +
  geom_point() +
  facet_wrap(. ~ measurement, ncol = 2)












#Is fout, wordt eigenlijk op de d_rvtlc data gedaan. Doen in for-loop

#Get TMM normalised counts
norm_counts <- cpm(d1, log = T)

#ENSG00000133110 ENSG00000161905 ENSG00000264940 
norm_counts %>% 
  as.data.frame() %>% 
  rownames_to_column("geneID") %>% 
  filter(geneID == "ENSG00000264940") %>% 
  pivot_longer(!geneID, names_to = "patient", values_to = "expr") %>% 
  mutate(patient = as.numeric(patient)) %>% 
  left_join(y = patient_data %>% select(idnr, d_ccq), by = c("patient" = "idnr")) %>% 
  ggplot(aes(y = expr, x = d_ccq)) + 
  geom_point() +
  labs(y = "Log2 TMM values", x = expression(Delta ~ "CCQ Total Score"))


for (i in top_genes){
  
  #Check for significant genes
  if (sum(i$FDR <= 0.05) > 0) {
    
    #Next, get TMM normalized counts
    norm_counts <- cpm(d1, log = T)
    
    #Select genes with FDR <= 0.05
    sign_genes <- i %>% 
      filter(FDR <= 0.05)
    
    for (gene in sign_genes$ensembl_gene_id) {
      print(gene)
      
      #plot
      norm_counts %>% 
        as.data.frame() %>% 
        rownames_to_column("geneID") %>% 
        filter(geneID == gene) %>% 
        pivot_longer(!geneID, names_to = "patient", values_to = "expr") %>% 
        mutate(patient = as.numeric(patient)) %>% 
        left_join(y = patient_data %>% select(idnr, d_ccq), by = c("patient" = "idnr"))
    }
    
    
    print("test")
    print(sign_genes)
  }
}
