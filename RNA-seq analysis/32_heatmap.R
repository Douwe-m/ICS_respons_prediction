library(RColorBrewer)
library(pheatmap)

GOI <- DEG_d_ccq %>% 
  filter(FDR <=0.10) %>% 
  pull(ensembl_gene_id)

heatmap_data <- norm_counts %>% 
  rownames_to_column("geneID") %>% 
  filter(geneID %in% GOI) %>% 
  left_join(gene_symbols %>% select(ensembl_gene_id, hgnc_symbol), by = c("geneID" = "ensembl_gene_id")) %>% 
  distinct() %>% 
  mutate(gene_symbol = ifelse(hgnc_symbol == "", geneID, hgnc_symbol)) %>% 
  select(-geneID, -hgnc_symbol) %>% 
  column_to_rownames("gene_symbol") %>% 
  as.matrix()

heatmap_annotation <- patient_data %>% 
  select(idnr, lftv1, treatment, rooknuv1) %>% 
  mutate(treatment = case_when(treatment == 2 ~ "ICS+LABA", 
                           treatment == 3 ~ "ICS",
                           treatment == 4 ~ "ICS+placebo")) %>% 
  mutate(rooknuv1 = ifelse(rooknuv1 == 0, "No", "Yes")) %>% 
  column_to_rownames("idnr")


pheatmap(heatmap_data,
         color = colorRampPalette(brewer.pal(8, "YlOrRd"))(25),
         annotation_col = heatmap_annotation)
