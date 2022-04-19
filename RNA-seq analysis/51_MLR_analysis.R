#Plot correlations
##DELTA FEV1 % pred
vars %>% 
  select(-d_rvtlc, -d_ccq) %>% 
  pivot_longer(!c(idnr, d_fev_pred), 
               names_to = "parameter", values_to = "measurement") %>% 
  ggplot(aes(x = measurement, y = d_fev_pred)) +
  geom_point() +
  facet_wrap(.~parameter, scales = "free")

##DELTA RV/TLC
vars %>% 
  select(-d_ccq, -d_fev_pred) %>% 
  pivot_longer(!c(idnr, d_rvtlc), 
               names_to = "parameter", values_to = "measurement") %>% 
  ggplot(aes(x = measurement, y = d_rvtlc)) +
  geom_point() +
  facet_wrap(.~parameter, scales = "free")

#DELTA CCQ total score
vars %>% 
  select(-d_fev_pred, -d_rvtlc) %>% 
  pivot_longer(!c(idnr, d_ccq), 
               names_to = "parameter", values_to = "measurement") %>% 
  ggplot(aes(x = measurement, y = d_ccq)) +
  geom_point() +
  facet_wrap(.~parameter, scales = "free")

#Calculate correlations
cor_matrix <-  corr.test(vars %>% select(-idnr), method = "spearman", adjust = "none")
cor_matrix_padj <-  corr.test(vars %>% select(-idnr), method = "spearman", adjust = "fdr")

correlations_r <- cor_matrix$r %>% as.data.frame() %>% select(starts_with("d_")) %>% rownames_to_column("parameter")
correlations_p <- cor_matrix$p %>% as.data.frame() %>% select(starts_with("d_")) %>% rownames_to_column("parameter")
correlations_padj <- cor_matrix_padj$p %>% as.data.frame() %>% select(starts_with("d_")) %>% rownames_to_column("parameter")

#Results from correlation analysis
d_fev_cor <- correlations_r %>% 
  select(parameter, d_fev_pred) %>% 
  filter(!grepl('d_', parameter)) %>% 
  left_join(correlations_p %>% select(parameter, d_fev_pred),
            by = "parameter") %>% 
  left_join(correlations_padj %>% select(parameter, d_fev_pred),
            by = "parameter") %>% 
  rename("r" = d_fev_pred.x, "p" = d_fev_pred.y, "FDR" = d_fev_pred) %>% 
  arrange(p)

d_ccq_cor <- correlations_r %>% 
  select(parameter, d_ccq) %>% 
  filter(!grepl('d_', parameter)) %>% 
  left_join(correlations_p %>% select(parameter, d_ccq),
            by = "parameter") %>% 
  left_join(correlations_padj %>% select(parameter, d_ccq),
            by = "parameter") %>% 
  rename("r" = d_ccq.x, "p" = d_ccq.y, "FDR" = d_ccq) %>% 
  arrange(p)

d_rvtlc_cor <- correlations_r %>% 
  select(parameter, d_rvtlc) %>% 
  filter(!grepl('d_', parameter)) %>% 
  left_join(correlations_p %>% select(parameter, d_rvtlc),
            by = "parameter") %>% 
  left_join(correlations_padj %>% select(parameter, d_rvtlc),
            by = "parameter") %>% 
  rename("r" = d_rvtlc.x, "p" = d_rvtlc.y, "FDR" = d_rvtlc) %>% 
  arrange(p) %>% 
  format(scientific=FALSE)

d_fev_cor
d_ccq_cor
d_rvtlc_cor

#multiple linear regression p<0.10
fit_d_ccq <- lm(formula = d_ccq ~ lftv1 + lymfopv1, data = patient_data_subset)
fit_d_fev <- lm(formula = d_fev_pred ~ lengtev1, data = patient_data_subset)
fit_d_rvtlc <- lm(formula = d_rvtlc ~ rvtlcpredv1 + fevfvcnapredv1 + ccqtotalv1 + monopv1 + lftv1 + fvcnapredv1 + 
                    eosv1, data = patient_data_subset)

#multiple linear regression p<0.20
fit_d_ccq <- lm(formula = d_ccq ~ lftv1 + lymfopv1 + ccqtotalv1 + signature + cysptlymv1 + neutropv1, data = patient_data_subset)
fit_d_fev <- lm(formula = d_fev_pred ~ lengtev1 + tlcv1 + cysptneuv1 + monov1, data = patient_data_subset)
fit_d_rvtlc <- lm(formula = d_rvtlc ~ rvtlcpredv1 + fevfvcnapredv1 + ccqtotalv1 + monopv1 + lftv1 + fvcnapredv1 + 
                    eosv1 + eospv1 + packyearv1 + spneutrov1 + cysptneuv1, data = patient_data_subset)

#multiple linear regression p<0.10 with genes
fit_d_ccq <- lm(formula = d_ccq ~ ENSG00000161905 + lftv1 + lymfopv1, data = vars)
fit_d_fev <- lm(formula = d_fev_pred ~ lengtev1 + ENSG00000172752, data = vars)
fit_d_rvtlc <- lm(formula = d_rvtlc ~ rvtlcpredv1 + fevfvcnapredv1 + ccqtotalv1 + monopv1 + lftv1 + fvcnapredv1 + 
                    eosv1, data = vars)

#multiple linear regression p<0.20 with genes
fit_d_ccq <- lm(formula = d_ccq ~ ENSG00000161905 + lftv1 + lymfopv1 + ccqtotalv1 + signature + cysptlymv1 + neutropv1, data = vars)
fit_d_fev <- lm(formula = d_fev_pred ~ lengtev1 + ENSG00000172752 + tlcv1 + cysptneuv1 + monov1, data = vars)
fit_d_rvtlc <- lm(formula = d_rvtlc ~ rvtlcpredv1 + fevfvcnapredv1 + ccqtotalv1 + monopv1 + lftv1 + fvcnapredv1 + 
                    eosv1 + eospv1 + packyearv1 + spneutrov1 + cysptneuv1, data = vars)

summary(fit_d_ccq)
summary(fit_d_fev)
summary(fit_d_rvtlc)

