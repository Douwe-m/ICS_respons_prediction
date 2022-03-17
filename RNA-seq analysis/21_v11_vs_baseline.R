#Filter samples available at baseline
counts_filtered <- counts %>% 
  select(ends_with("_1")) %>% 
  rename_with(~gsub("_1|X", "", .x))

#Select required columns
patient_data <- patient_data %>% 
  select(idnr, treatment, geslacht, lftv1, rookjaar, rooknuv1, packyearv1, 
         ccqtotalv1, ccqtotalv11, rvtlcpredv1, rvtlcpredv11, fevnapredv1, fevnapredv11) 

#Filter for patients that have RNA-seq available at baseline, no NA in any column and have received no placebo
patient_data <- patient_data %>% 
  filter(idnr %in% as.numeric(colnames(counts_filtered))) %>% 
  filter(treatment != 1) %>%
  filter(treatment != 4) %>%
  drop_na()

#Add a column with the change in lung function after 30 months
patient_data <- patient_data %>% 
  mutate(d_fev_pred = fevnapredv11 - fevnapredv1) %>% 
  mutate(d_ccq = ccqtotalv11 - ccqtotalv1) %>% 
  mutate(d_rvtlc = rvtlcpredv11 - rvtlcpredv1)

#Final selection of patients
counts_filtered <- counts_filtered %>%
  select(which(colnames(counts_filtered) %in% patient_data$idnr))

#Create an empty list to add results to
res_v11 <- list()

#Create a DGEList object
d0 <- DGEList(counts_filtered)

#Filter lowly expressed genes
keep <- filterByExpr(d0)
d0 <- d0[keep, , keep.lib.sizes = F]

#Normalize counts
d1 <- calcNormFactors(d0)

#Get TMM normalized counts
norm_counts <- cpm(d1, log = T)

clin_data <- patient_data %>% select(idnr, starts_with("d_"))

for (i in 2:ncol(clin_data)){
  print(colnames(clin_data[,i]))
  clin_param <- clin_data %>% pull(i)
  
  #Design matrix
  mm <- model.matrix(~clin_param)
  row.names(mm) <- clin_data$idnr
  
  #estimate dispersion
  d2 <- estimateDisp(d1, mm, robust = T)
  
  fit <- glmQLFit(d2, mm)
  fit <- glmQLFTest(fit)
  
  #Get results
  tt <- topTags(fit, n = nrow(d2), adjust.method = "BH")
  
  #Add ensemble gene IDs
  top_tt <- tt$table %>%
    rownames_to_column("ensembl_gene_id") %>%
    left_join(y = gene_symbols %>% select(ensembl_gene_id, hgnc_symbol), by = "ensembl_gene_id") %>% 
    distinct()
  
  #Add findings to the named list of results
  res_v11 <- c(list(top_tt), res_v11)
}











