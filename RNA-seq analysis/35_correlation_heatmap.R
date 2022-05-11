library(ggcorrplot)

GOI <- DEG_d_ccq %>% 
  filter(FDR <= 0.10) %>% 
  pull(ensembl_gene_id)

DEG_norm_counts <- norm_counts %>% 
  filter(rownames(.) %in% GOI) %>% 
  rownames_to_column("geneID") %>% 
  left_join(gene_symbols %>% select(ensembl_gene_id, hgnc_symbol), by = c("geneID" = "ensembl_gene_id")) %>% 
  mutate(geneID = ifelse(hgnc_symbol == "" | is.na(hgnc_symbol), 
                         geneID,
                         hgnc_symbol)) %>% 
  select(-hgnc_symbol) %>% 
  distinct() %>% 
  column_to_rownames("geneID") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("patientID") %>% 
  mutate(patientID = as.numeric(patientID))

es <- gsva_es %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("patientID") %>%
  mutate(patientID = as.numeric(patientID))

patient_data_subset <- patient_data_raw %>% 
  filter(treatment != 1) %>%
  left_join(es, by = c("idnr" = "patientID")) %>% 
  select(idnr, treatment, lftv1, packyearv1, lengtev1, gewichtv1,   
         ccqtotalv1, ccqtotalv3, ccqtotalv11,
         rvtlcpredv1, rvtlcpredv3, rvtlcpredv11,
         fevnapredv1, fevnapredv3, fevnapredv11,
         fvcnapredv1, fevfvcnapredv1, tlcv1,
         eospv1, eosv1, speov1, cyspteov1, 
         neutropv1, spneutrov1, cysptneuv1,
         lymfopv1, lymfov1, splymv1, cysptlymv1,
         monopv1, monov1, 
         signature = d_ccq_sig) %>% 
  left_join(DEG_norm_counts, by = c("idnr" = "patientID")) %>% 
  mutate(d_fev_pred = fevnapredv3 - fevnapredv1,
         d_ccq = ccqtotalv3 - ccqtotalv1,
         d_rvtlc = rvtlcpredv3 - rvtlcpredv1,
         v11_d_fev_pred = fevnapredv11 - fevnapredv1,
         v11_d_ccq = ccqtotalv11 - ccqtotalv1,
         v11_d_rvtlc = rvtlcpredv11 - rvtlcpredv1)

corr_data_6_months <- patient_data_subset %>% 
  select(-ends_with("v11"), -ends_with("v3"), -starts_with("v11_"), -treatment, -idnr) %>% 
  drop_na()

corr_data_30_months <- patient_data_subset %>% 
  filter(treatment != 4) %>% 
  select(-ends_with("v11"), -ends_with("v3"), -starts_with("d_"), -treatment, -idnr) %>% 
  drop_na()


cor_matrix_v1 <-  corr.test(corr_data_6_months, method = "spearman", adjust = "fdr")
correlations_r_v1 <- cor_matrix_v1$r %>% as.data.frame() %>% select(starts_with("d_")) %>% rownames_to_column("parameter")
correlations_padj_v1 <- cor_matrix_v1$p %>% as.data.frame() %>% select(starts_with("d_")) %>% rownames_to_column("parameter")

cor_matrix_v11 <-  corr.test(corr_data_30_months, method = "spearman", adjust = "fdr")
correlations_r_v11 <- cor_matrix_v11$r %>% as.data.frame() %>% select(starts_with("v11_")) %>% rownames_to_column("parameter")
correlations_padj_v11 <- cor_matrix_v11$p %>% as.data.frame() %>% select(starts_with("v11_")) %>% rownames_to_column("parameter")

correlations_r <- correlations_r_v1 %>% 
  full_join(correlations_r_v11, by = "parameter") %>% 
  filter(str_detect(parameter, "d_", negate = T)) %>% 
  column_to_rownames("parameter")

correlations_padj <- correlations_padj_v1 %>% 
  full_join(correlations_padj_v11, by = "parameter") %>% 
  filter(str_detect(parameter, "d_", negate = T)) %>% 
  column_to_rownames("parameter")

ggcorrplot(correlations_r, colors = c("#6D9EC1", "white", "#E46726"))
