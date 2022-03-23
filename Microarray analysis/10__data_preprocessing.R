#Add gene names to expression data
expression_data_cleaned <- expression_data %>% 
  left_join(gene_names %>% 
              select(rowid, SYMBOL),
            by = "rowid") %>% 
  select(-rowid, symbol = SYMBOL) %>% 
  relocate(symbol)

#Replace column names with patient IDs
#Thanks to: https://stackoverflow.com/questions/60113369/change-column-names-in-dataframe-based-on-matching-to-another-dataframe-by-dplyr
IDs <- patient_data %>% 
  select(1, 2, 5) %>% 
  mutate(id = paste0(patientID, "_V", time)) %>% 
  select(id, rowid) %>% 
  deframe()

expression_data_cleaned <- expression_data_cleaned %>% 
  rename(!!!IDs) %>% 
  replace_na(list(symbol = "unknown_gene")) %>% #Replace the one NA in the symbol column with "unknown_gene"
  column_to_rownames("symbol")

#Filter expression data at baseline
expression_data_cleaned <- expression_data_cleaned %>% 
  select(ends_with("_V0")) %>% 
  rename_with(~gsub("_V0", "", .x))

#Clean patient data
patient_data <- patient_data %>% 
  select(patientID, treatment, time, gender, agev1, smokingstatusv1, packyearv1, 
         ccqtotalv1, ccqtotalv3, rvtlcpredv1, rvtlcpredv3, fevnaprev1, fevnaprev3) %>% 
  filter(treatment != 0,
         time == 0) %>% 
  mutate(d_fev_pred = fevnaprev3 - fevnaprev1,
         d_ccq = ccqtotalv3 - ccqtotalv1,
         d_rvtlc = rvtlcpredv3 - rvtlcpredv1) %>% 
  drop_na()

#Final selection of patients
expression_data_cleaned <- expression_data_cleaned %>%
  select(which(colnames(expression_data_cleaned) %in% patient_data$patientID))


