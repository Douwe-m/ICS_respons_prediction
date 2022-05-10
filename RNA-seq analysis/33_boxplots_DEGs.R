boxplot_DEGs <- ggplot(patient_data_visualisation, aes(x = geneID, y = cpm_counts)) +
  geom_jitter(width = 0.25, height = 0) +
  geom_boxplot(fill = "lightblue", outlier.shape = NA) +
  scale_y_continuous(breaks = seq(10, 20, 1), limits = c(10, 20)) +
  labs(x = "", y = "log<sub>2</sub>(CPM)") +
  my_theme +
  theme(axis.title.x = element_markdown()) +
  coord_flip()
    
boxplot_DEGs
