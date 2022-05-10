#Theme for data visualizations 
my_theme<-theme_classic()+
  theme(panel.border = element_rect(color='black', fill=NA, size=1.2),
        axis.ticks = element_line(color='dimgray', size=0.8),
        axis.text = element_text(size=9,colour = 'dimgray'),
        axis.title = element_text(size=14),
        plot.title = element_text(size=17),
        legend.title = element_text(size=14), 
        legend.position = 'bottom', 
        text = element_text(family = "serif"))

#Select significant DEGs
GOI <- DEG_d_ccq %>% 
  filter(FDR <= 0.10) %>% 
  pull(ensembl_gene_id)

#Add proper gene IDs
toptags_d_ccq <- DEG_d_ccq %>% 
  mutate(geneID = ifelse(hgnc_symbol == "" | is.na(hgnc_symbol), 
                         ensembl_gene_id,
                         hgnc_symbol)) %>% 
  select(-ensembl_gene_id, -hgnc_symbol)

#CPM counts of DEGs
cpm_counts_DEGs <- counts_filtered %>% 
  filter(rownames(.) %in% GOI) %>% 
  cpm(log = T) %>% 
  as.data.frame() %>% 
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

patient_data_visualisation <- patient_data %>% 
  left_join(cpm_counts_DEGs, by = c("idnr" = "patientID")) %>% 
  select(idnr, contains("d_"), contains("ENSG")) %>% 
  pivot_longer(contains("ENSG"), names_to = "geneID", values_to = "cpm_counts") %>% 
  left_join(gene_symbols %>% select(ensembl_gene_id, hgnc_symbol), by = c("geneID" = "ensembl_gene_id")) %>% 
  mutate(geneID = ifelse(hgnc_symbol == "", geneID, hgnc_symbol)) %>% 
  select(-hgnc_symbol) %>% 
  distinct()

