p1 <- volcano_plot
p2 <- as.ggplot(heatmap)
p3 <- boxplot_DEGs


(p1+p3)/p2 + plot_annotation(tag_levels = "A")
