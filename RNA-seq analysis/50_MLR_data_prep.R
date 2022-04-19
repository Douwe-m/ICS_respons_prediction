#prepare enrichment scores
es <- gsva_es %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("patient") %>% 
  mutate(patient = as.numeric(patient))

#Get CPM norm counts
GOI <- DEG_d_ccq %>% 
  filter(FDR <= 0.10) %>% 
  pull(ensembl_gene_id)

cpm_counts <- cpm(counts_filtered, log = T) %>% 
  as.data.frame() %>% 
  rownames_to_column("geneID") %>% 
  filter(geneID %in% GOI) %>% 
  pivot_longer(!geneID, names_to = "patient", values_to = "expr") %>% 
  pivot_wider(names_from = geneID, values_from = expr) %>% 
  mutate(patient = as.numeric(patient))

#clean patient data
patient_data_subset <- patient_data_raw %>% 
  filter(treatment != 1) %>% 
  left_join(es, by = c("idnr" = "patient")) %>% 
  select(idnr, treatment, geslacht, lftv1, rooknuv1, packyearv1, lengtev1, gewichtv1,   
         ccqtotalv1, ccqtotalv3, rvtlcpredv1, rvtlcpredv3, fevnapredv1, fevnapredv3,
         fvcnapredv1, fevfvcnapredv1, tlcv1,
         eospv1, eosv1, speov1, cyspteov1, 
         neutropv1, spneutrov1, cysptneuv1,
         lymfopv1, lymfov1, splymv1, cysptlymv1,
         monopv1, monov1, 
         signature = d_ccq_sig) %>% 
  mutate(d_fev_pred = fevnapredv3 - fevnapredv1,
         d_ccq = ccqtotalv3 - ccqtotalv1,
         d_rvtlc = rvtlcpredv3 - rvtlcpredv1) %>% 
  left_join(cpm_counts, by = c("idnr" = "patient"))

#Remove variables not to include
vars <- patient_data_subset %>% 
  drop_na() %>% 
  select(-treatment, -geslacht, -rooknuv1, -ccqtotalv3, -rvtlcpredv3, -fevnapredv3)
