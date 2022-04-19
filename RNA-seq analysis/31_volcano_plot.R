library(EnhancedVolcano)
library(ggrepel)

my_theme<-theme_classic()+
  theme(panel.border = element_rect(color='black', fill=NA, size=1.2),
        axis.ticks = element_line(color='dimgray', size=0.8),
        axis.text = element_text(size=9,colour = 'dimgray'),
        axis.title = element_text(size=14),
        plot.title = element_text(size=17),
        legend.title = element_text(size=14), 
        legend.position = 'bottom')

EnhancedVolcano(DEG_d_ccq, 
                x = "logFC",
                y = "FDR",
                lab = "ensembl_gene_id",
                pCutoff = 0.10)


new_DEGs <- DEG_d_ccq %>% 
  mutate(temp_symbols = ifelse(hgnc_symbol == "", NA, hgnc_symbol)) %>% 
  mutate(gene_symbol = ifelse(is.na(temp_symbols), ensembl_gene_id, temp_symbols)) %>% 
  select(-hgnc_symbol, -temp_symbols)


ggplot(new_DEGs, aes(x = logFC, y = -log10(FDR))) +
  geom_hline(aes(yintercept = -log10(0.10)), linetype = "dashed") +
  geom_point(alpha = 0.25) +
  geom_point(data = subset(new_DEGs, FDR<=0.10), color = "lightcoral", size = 3) +
  geom_text_repel(data = subset(new_DEGs, FDR<=0.10), aes(x = logFC, y = -log10(FDR), label = gene_symbol)) +
  scale_x_continuous(breaks = seq(-3, 3, 1), limits = c(-3, 3)) +
  scale_y_continuous(breaks = seq(0, 3, 1), limits = c(0, 3)) +
  my_theme
