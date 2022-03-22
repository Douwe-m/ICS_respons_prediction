#Scatter plot of expression vs change in clinical parameter

if (sum(DEG_d_ccq$FDR <= 0.25) > 0) {

  scatter_plot <- norm_counts %>% 
    rename_with(~gsub("X", "", .x)) %>% 
    rownames_to_column("geneID") %>% 
    filter(geneID %in% DEG_d_ccq[DEG_d_ccq$FDR <=0.05, 1]) %>% 
    pivot_longer(!geneID, names_to = "patient", values_to = "expr") %>% 
    mutate(patient = as.numeric(patient)) %>%
    left_join(y = patient_data %>% select(idnr, d_ccq), by = c("patient" = "idnr")) %>% 
    
    ggplot(aes(x = expr, y = d_ccq, group = factor(geneID))) +
      geom_point() +
      labs(x = "Log2 TMM values", y = expression(Delta ~ "CCQ Total Score")) +
      facet_wrap(. ~ geneID, ncol = 2) 
  
  ggsave(filename = "RNA-seq analysis/results/scatter_plot_d_ccq.png",
         plot = scatter_plot)
}

if (sum(DEG_d_fev$FDR <= 0.25) > 0) {
  
  scatter_plot <- norm_counts %>% 
    rename_with(~gsub("X", "", .x)) %>% 
    rownames_to_column("geneID") %>% 
    filter(geneID %in% DEG_d_fev[DEG_d_fev$FDR <=0.05, 1]) %>% 
    pivot_longer(!geneID, names_to = "patient", values_to = "expr") %>% 
    mutate(patient = as.numeric(patient)) %>%
    left_join(y = patient_data %>% select(idnr, d_fev_pred), by = c("patient" = "idnr")) %>% 
    
    ggplot(aes(x = expr, y = d_fev_pred, group = factor(geneID))) +
      geom_point() +
      labs(x = "Log2 TMM values", y = expression(Delta ~ "FEV1 % pred")) +
      facet_wrap(. ~ geneID, ncol = 2) 
  
  ggsave(filename = "RNA-seq analysis/results/scatter_plot_fev.png",
         plot = scatter_plot)
}

if (sum(DEG_d_rvtlc$FDR <= 0.25) > 0) {
  
  scatter_plot <- norm_counts %>% 
    rename_with(~gsub("X", "", .x)) %>% 
    rownames_to_column("geneID") %>% 
    filter(geneID %in% DEG_d_rvtlc[DEG_d_rvtlc$FDR <=0.05, 1]) %>% 
    pivot_longer(!geneID, names_to = "patient", values_to = "expr") %>% 
    mutate(patient = as.numeric(patient)) %>%
    left_join(y = patient_data %>% select(idnr, d_rvtlc), by = c("patient" = "idnr")) %>% 
    
    ggplot(aes(x = expr, y = d_rvtlc, group = factor(geneID))) +
      geom_point() +
      labs(x = "Log2 TMM values", y = expression(Delta ~ "RV/TLC % pred")) +
      facet_wrap(. ~ geneID, ncol = 2) 
  
  ggsave(filename = "RNA-seq analysis/results/scatter_plot_d_rvtlc.png",
         plot = scatter_plot)
}

#3x3 scatter plot of genes of interest (GOI)
GOI <- c("ENSG00000133110", "ENSG00000161905", "ENSG00000264940")

temp_data <- norm_counts %>% 
  rename_with(~gsub("X", "", .x)) %>% 
  rownames_to_column("geneID") %>% 
  filter(geneID %in% GOI) %>% 
  pivot_longer(!geneID, names_to = "patient", values_to = "expr") %>% 
  mutate(patient = as.numeric(patient)) %>%
  left_join(y = patient_data %>% select(idnr, d_ccq, d_fev_pred, d_rvtlc), by = c("patient" = "idnr")) %>% 
  pivot_longer(cols = c(d_ccq, d_fev_pred, d_rvtlc), names_to = "clin_param", values_to = "values")

plot <- temp_data %>%   
  ggplot(aes(x = expr, y = values)) +
    geom_point() +
    geom_smooth(method = lm, se = F)+
    facet_wrap(geneID ~ clin_param, scales = "free") +
    geom_label(aes(label = patient))
  
ggsave(filename = "RNA-seq analysis/results/scatter_plot_3_genes.png",
       plot = plot)

#Calculate correlations
temp_data %>% 
  group_by(geneID, clin_param) %>% 
  summarise(COR = cor.test(expr, values, method = "spearman", exact = F)$estimate, 
            pval = cor.test(expr, values, method = "spearman", exact = F)$p.value) %>% 
  print()









