#Plot change in clinical parameters with enrichment score after 6 months
signature_corr_6months <- ggplot(gsva_es_cleaned, aes(x = enrichment_score, y = value)) +
  geom_point() +
  # geom_smooth(method = "lm", se = F) +
  stat_cor(method = "spearman") +
  labs(x = "Enrichment score", y = "Clinical parameter") +
  scale_x_continuous(breaks = seq(-0.5, 1, 0.5), limits = c(-0.75, 1)) +
  my_theme +
  facet_wrap(. ~ parameter, ncol = 3, scales = "free",
             labeller = labeller(parameter = 
                                   c("d_ccq" = "ΔCCQ",
                                     "d_fev_pred" = "ΔFEV1 % pred.",
                                     "d_rvtlc" = "ΔRV/TLC % pred."))) +
  theme(aspect.ratio = 1)
  

#Plot change in clinical parameters with enrichment score after 30 months
signature_corr_30months <- ggplot(gsva_es_cleaned_v11, aes(x = enrichment_score, y = value)) +
  geom_point() +
  # geom_smooth(method = "lm", se = F) +
  stat_cor(method = "spearman") +
  labs(x = "Enrichment score", y = "Clinical parameter") +
  scale_x_continuous(breaks = seq(-0.5, 1, 0.5), limits = c(-0.75, 0.6)) +
  my_theme +
  facet_wrap(. ~ parameter, ncol = 3, scales = "free",
             labeller = labeller(parameter = 
                                   c("d_ccq" = "ΔCCQ",
                                     "d_fev_pred" = "ΔFEV1 % pred.",
                                     "d_rvtlc" = "ΔRV/TLC % pred."))) +
  theme(aspect.ratio = 1)

