library(EnhancedVolcano)
library(ggrepel)
library("ggtext")

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

y.bar.intercept <- new_DEGs %>%
  dplyr::filter(FDR >=0.10) %>% 
  dplyr::pull(PValue) %>%
  min() %>%
  log10() %>%
  (function(x) {x * -1})



ggplot(new_DEGs, aes(x = logFC, y = -log10(PValue))) +
  geom_hline(aes(yintercept = y.bar.intercept), linetype = "dashed") +
  geom_point(alpha = 0.25) +
  geom_point(data = subset(new_DEGs, FDR<=0.10), color = "lightcoral", size = 3) +
  geom_text_repel(data = subset(new_DEGs, FDR<=0.10), aes(x = logFC, y = -log10(PValue), label = gene_symbol)) +
  scale_x_continuous(breaks = seq(-3, 3, 1), limits = c(-3, 3)) +
  scale_y_continuous(breaks = seq(0, 8, 1), limits = c(0, 8)) +
  my_theme +
  theme(axis.title.y = element_markdown(),
        axis.title.x = element_markdown()) +
  labs(x = "log<sub>2</sub>(FC)", 
       y = "-log<sub>10</sub>(p)")






