# if (sum(DEG_d_fev$FDR <= 0.25) > 0){}
DEG_d_ccq %>% 
  filter(adj.P.Val <= 0.05) %>% 
  count()

DEG_d_fev %>% 
  filter(adj.P.Val <= 0.05) %>% 
  count()

DEG_d_rvtlc %>% 
  filter(adj.P.Val <= 0.05) %>% 
  count()

#Genes of interest
gene_set <- list(
  d_ccq_all = DEG_d_ccq %>% 
    filter(adj.P.Val <= 0.05) %>% 
    rownames_to_column("geneID") %>% 
    pull(geneID),
  d_ccq_up = DEG_d_ccq %>% 
    filter(adj.P.Val <= 0.05,
           logFC > 0) %>% 
    rownames_to_column("geneID") %>% 
    pull(geneID),
  d_ccq_down = DEG_d_ccq %>% 
    filter(adj.P.Val <= 0.05,
           logFC < 0) %>% 
    rownames_to_column("geneID") %>% 
    pull(geneID))

data <- voom(expression_data_cleaned)

#Calculate GSVA enrichment scores
gsva_es <- gsva(as.matrix(data), gene_set)

gsva_es_cleaned <- gsva_es %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("patient") %>% 
  mutate(patient = as.numeric(patient)) %>% 
  left_join(patient_data %>% 
              select(patientID, d_ccq, d_fev_pred, d_rvtlc),
            by = c("patient" = "patientID")) %>% 
  pivot_longer(c(d_ccq_all, d_ccq_up, d_ccq_down), names_to = "signature", values_to = "enrichment_score") %>% 
  pivot_longer(c(d_ccq, d_fev_pred, d_rvtlc), names_to = "clin_param", values_to = "value")

#Plot change in clinical parameters with enrichment score
gsva_es_cleaned %>% 
  ggplot(aes(x = enrichment_score, y = value)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Enrichment score", y = "Clinical parameter") +
  facet_wrap(signature ~ clin_param, scales = "free") 

#Calculate correlations
gsva_es_cleaned %>% 
  group_by(signature, clin_param) %>% 
  summarise(COR = cor.test(enrichment_score , value, method = "spearman", exact = F)$estimate, 
            pval = cor.test(enrichment_score, value, method = "spearman", exact = F)$p.value) %>% 
  print()
