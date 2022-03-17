#Remove patients 2109 and 2107 since they are possible outliers 
counts_pat_rem <- counts_filtered %>% 
  select(-`2109`, -`2107`)

patient_data_pat_rem <- patient_data %>% 
  filter(idnr != 2109,  
         idnr != 2107)

#Create an empty list to add results to
top_genes <- list()

#Create a DGEList object
d0 <- DGEList(counts_pat_rem)

#Filter lowly expressed genes
keep <- filterByExpr(d0)
d0 <- d0[keep, , keep.lib.sizes = F]

#Normalize counts
d1 <- calcNormFactors(d0)

#Get TMM normalized counts
norm_counts <- cpm(d1, log = T)

clin_data <- patient_data_pat_rem %>% select(idnr, starts_with("d_"))

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
  top_genes <- c(list(top_tt), top_genes)
}


