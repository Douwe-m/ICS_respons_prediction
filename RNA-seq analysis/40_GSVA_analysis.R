#GSVA requires two things:
#1) A normalized gene expression dataset
#2) A collection of gene sets as a list

#Genes of interest
gene_set <- list("d_ccq_sig" = DEG_d_ccq %>% 
                   filter(FDR <= 0.10) %>% 
                   pull(ensembl_gene_id))

#Calculate GSVA enrichment scores
gsva_es <- gsva(as.matrix(norm_counts), gene_set)

gsva_es_cleaned <- gsva_es %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("patient") %>%
  rename(enrichment_score = d_ccq_sig) %>%
  mutate(patient = as.numeric(patient)) %>%
  left_join(patient_data %>%
              select(idnr, d_ccq, d_fev_pred, d_rvtlc),
            by = c("patient" = "idnr")) %>%
  pivot_longer(3:5, names_to = "parameter", values_to = "value")

#Calculate correlations
gsva_es_cleaned %>%
  group_by(parameter) %>%
  summarise(COR = cor.test(enrichment_score , value, method = "spearman", exact = F)$estimate,
            pval = cor.test(enrichment_score, value, method = "spearman", exact = F)$p.value) %>%
  print()


#Calculate GSVA enrichment scores
gsva_es_v11 <- gsva(as.matrix(norm_counts_v11), gene_set)

gsva_es_cleaned_v11 <- gsva_es_v11 %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("patient") %>%
  rename(enrichment_score = d_ccq_sig) %>%
  mutate(patient = as.numeric(patient)) %>%
  left_join(patient_data_v11 %>%
              select(idnr, d_ccq, d_fev_pred, d_rvtlc),
            by = c("patient" = "idnr")) %>%
  pivot_longer(3:5, names_to = "parameter", values_to = "value") %>%
  drop_na()

#Calculate correlations
gsva_es_cleaned_v11 %>%
  group_by(parameter) %>%
  summarise(COR = cor.test(enrichment_score , value, method = "spearman", exact = F)$estimate,
            pval = cor.test(enrichment_score, value, method = "spearman", exact = F)$p.value) %>%
  print()

