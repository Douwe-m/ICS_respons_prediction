GOI <- DEG_d_ccq %>% 
  filter(FDR <= 0.10) %>% 
  pull(ensembl_gene_id)

#Get CPM norm counts
cpm_counts <- cpm(counts_filtered, log = T) %>% 
  as.data.frame() %>% 
  rownames_to_column("geneID") %>% 
  filter(geneID %in% GOI) %>% 
  pivot_longer(!geneID, names_to = "patient", values_to = "expr") %>% 
  pivot_wider(names_from = geneID, values_from = expr) %>% 
  mutate(patient = as.numeric(patient))
  
#Add CPM for genes of interest
patient_data_subset <- patient_data_subset %>% 
  left_join(cpm_counts, by = c("idnr" = "patient"))

#Drop columns that will not be used in the analysis
vars <- patient_data_subset %>% 
  select(-treatment, -geslacht, -rookjaar, -rooknuv1, -ccqtotalv3, -rvtlcpredv3, -fevnapredv3) %>% 
  drop_na()
  