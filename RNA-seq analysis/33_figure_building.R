#Figure 1
p1 <- volcano_plot
p2 <- as.ggplot(heatmap)
p3 <- boxplot_DEGs

(p1+p3)/p2 + plot_annotation(tag_levels = "A")

#Figure 2
p4 <- signature_corr_6months
p5 <- signature_corr_30months

p4/p5 + plot_annotation(tag_levels = "A")
