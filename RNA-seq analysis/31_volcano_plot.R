#Determine the intercept of the p-value cutoff
y_bar_intercept <- toptags_d_ccq %>%
  filter(FDR >=0.10) %>% 
  pull(PValue) %>%
  min() %>%
  log10() %>%
  (function(x) {x * -1}) 

volcano_plot <- ggplot(toptags_d_ccq, aes(x = logFC, y = -log10(PValue))) +
  geom_hline(aes(yintercept = y_bar_intercept), linetype = "dashed") +
  geom_point(alpha = 0.25) +
  geom_point(data = subset(toptags_d_ccq, FDR<=0.10), 
             color = "lightcoral", size = 3) +
  geom_text_repel(data = subset(toptags_d_ccq, FDR<=0.10), 
                  aes(x = logFC, y = -log10(PValue), label = geneID),
                  family = "serif") +
  scale_x_continuous(breaks = seq(-3, 2, 1), limits = c(-3, 2.1)) +
  scale_y_continuous(breaks = seq(0, 8, 1), limits = c(0, 8)) +
  my_theme +
  theme(axis.title.y = element_markdown(),
        axis.title.x = element_markdown()) +
  labs(x = "log<sub>2</sub>(FC)", 
       y = "-log<sub>10</sub>(p)")

volcano_plot
