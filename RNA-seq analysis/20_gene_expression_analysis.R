#Create a DGEList object
d0 <- DGEList(counts_filtered)

#Filter lowly expressed genes
keep <- filterByExpr(d0)
d0 <- d0[keep, , keep.lib.sizes = F]

#Normalize counts
d1 <- calcNormFactors(d0)

#Get TMM normalized counts
norm_counts <- cpm(d1, log = T)
write.table(norm_counts, file = "RNA-seq analysis/output/normalised_counts.csv", sep = ",")

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
  
  write.table(x = top_tt, 
              file = paste0("RNA-seq analysis/output/DEG_", colnames(clin_data[,i]), "_v3_baseline", ".csv"), 
              sep = ",")

}
