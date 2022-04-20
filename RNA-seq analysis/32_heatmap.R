library(RColorBrewer)

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

heatmap(heatmap_data, col = colorRampPalette(brewer.pal(8, "YlOrRd"))(25))
