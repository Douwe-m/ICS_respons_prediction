#Correlate identified genes to change in clinical parameters
GOI <- DEG_d_ccq %>% 
  filter(FDR <= 0.10) %>% 
  pull(ensembl_gene_id)

temp_data <- norm_counts %>% 
  rownames_to_column("geneID") %>% 
  filter(geneID %in% GOI) %>% 
  pivot_longer(!geneID, names_to = "patient", values_to = "expr") %>% 
  mutate(patient = as.numeric(patient)) %>%
  left_join(y = patient_data %>% select(idnr, d_ccq, d_fev_pred, d_rvtlc), by = c("patient" = "idnr")) %>% 
  pivot_longer(cols = c(d_ccq, d_fev_pred, d_rvtlc), names_to = "clin_param", values_to = "values")

temp_data %>%   
  ggplot(aes(x = expr, y = values)) +
  geom_point() +
  geom_smooth(method = lm, se = F)+
  facet_wrap(clin_param ~ geneID, scales = "free", ncol = 6)

# ggsave(filename = "RNA-seq analysis/results/scatter_plot_3_genes.png",
#        plot = plot)

#Calculate correlations
temp_data %>% 
  group_by(geneID, clin_param) %>% 
  summarise(COR = cor.test(expr, values, method = "spearman", exact = F)$estimate, 
            pval = cor.test(expr, values, method = "spearman", exact = F)$p.value) %>% 
  print()



temp_data %>% 
  select(geneID, patient, expr) %>% 
  distinct() %>% 
  ggplot(aes(x = geneID, y = expr)) +
  geom_boxplot() +
  geom_jitter()
   

temp_data %>% 
  select(geneID, patient, expr) %>% 
  distinct() %>% 
  ggplot(aes(x = geneID, y = expr)) +
  geom_boxplot() +
  geom_jitter()