#Function to calculate z-scores
calc_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

#Data to make heatmap with
heatmap_data <- cpm_counts_DEGs %>% 
  column_to_rownames("patientID") %>% 
  apply(2, calc_z_score) %>% 
  t()

#Heatmap annotation data
heatmap_annotation <- patient_data %>% 
  select(idnr, "Treatment" = treatment, "Smoking status" = rooknuv1, "ΔCCQ" = d_ccq) %>% 
  mutate(Treatment = case_when(Treatment == 2 ~ "ICS+LABA", 
                               Treatment == 3 ~ "ICS",
                               Treatment == 4 ~ "ICS+placebo")) %>% 
  mutate(`Smoking status` = ifelse(`Smoking status` == 0, "Former", "Current")) %>% 
  column_to_rownames("idnr") %>% 
  arrange(ΔCCQ)

#Color scale for the ΔCCQ annotation
annotation_colors = list(ΔCCQ = colorRampPalette(brewer.pal(6, "Blues"))(25))

#Draw the heatmap
heatmap <- pheatmap(heatmap_data[, rownames(heatmap_annotation)],
                     color = colorRampPalette(brewer.pal(8, "YlOrRd"))(25),
                     annotation_col = heatmap_annotation, 
                     cluster_cols = F,
                     cluster_rows = F, 
                     annotation_colors = annotation_colors,
                     show_colnames = F)

