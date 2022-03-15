#Get genes with an FDR<0.05
sign_genes <- data.frame()

for (i in top_genes) {
  if (sum(i$FDR <= 0.05) > 0){
    sign_genes <- rbind(sign_genes, i %>% filter(FDR<=0.05))
  }
}

#create df with counts of sign. genes and clinical parameter
a <- norm_counts %>% 
  as.data.frame() %>% 
  rownames_to_column("geneID") %>% 
  filter(geneID %in% sign_genes$ensembl_gene_id) %>% 
  pivot_longer(2:37, names_to = "patient", values_to = "expr") %>% 
  mutate(patient = as.numeric(patient)) %>%
  left_join(y = patient_data %>% select(idnr, d_ccq), by = c("patient" = "idnr"))

#scatter plot
ggplot(a, aes(y = d_ccq, x = expr, group = factor(geneID))) + 
  geom_point() +
  labs(x = "Log2 TMM values", y = expression(Delta ~ "CCQ Total Score"))+
  facet_wrap(. ~ geneID, ncol = 2)

#3x3
a %>%
  left_join(patient_data %>% select(idnr, d_fev_pred, d_rvtlc), by = c("patient" = "idnr")) %>% 
  pivot_longer(cols = c(d_ccq, d_fev_pred, d_rvtlc), names_to = "clin_param", values_to = "values") %>% 
  ggplot(aes(x = expr, y = values)) +
  geom_point() +
  geom_smooth(method = lm, se = F)+
  facet_wrap(geneID ~ clin_param, scales = "free")
