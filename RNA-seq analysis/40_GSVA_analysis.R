#GSVA requires two things:
#1) A normalized gene expression dataset
#2) A collection of gene sets as a list

#Genes of interest
gene_set <- list("d_ccq" = c("ENSG00000133110", "ENSG00000161905", "ENSG00000264940"))

#Calculate GSVA enrichment scores
gsva_es <- gsva(as.matrix(norm_counts), gene_set)

gsva_es_cleaned <- gsva_es %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("patient") %>% 
  rename(enrichment_score = d_ccq) %>% 
  mutate(patient = as.numeric(patient)) %>% 
  left_join(patient_data %>% 
              select(idnr, d_ccq, d_fev_pred, d_rvtlc),
            by = c("patient" = "idnr")) %>% 
  pivot_longer(3:5, names_to = "parameter", values_to = "value")

#Plot change in clinical parameters with enrichment score
gsva_es_cleaned %>% 
  ggplot(aes(x = enrichment_score, y = value, group = factor(parameter))) +
    geom_point() +
    geom_smooth(method = "lm") +
    labs(x = "Enrichment score", y = "Clinical parameter") +
    facet_wrap(. ~ parameter, ncol = 2, scales = "free") 

#Calculate correlations
gsva_es_cleaned %>% 
  group_by(parameter) %>% 
  summarise(COR = cor.test(enrichment_score , value, method = "spearman", exact = F)$estimate, 
            pval = cor.test(enrichment_score, value, method = "spearman", exact = F)$p.value) %>% 
  print()
