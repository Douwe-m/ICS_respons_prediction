#Filter samples available at baseline
counts_filtered <- counts %>% 
  select(ends_with("_1")) %>% 
  rename_with(~gsub("_1|X", "", .x))

#Select required columns
patient_data <- patient_data_raw %>% 
  select(idnr, treatment, geslacht, lftv1, rookjaar, rooknuv1, packyearv1, 
         ccqtotalv1, ccqtotalv3, rvtlcpredv1, rvtlcpredv3, fevnapredv1, fevnapredv3) 

#Filter for patients that have RNA-seq available at baseline, no NA in any column and have received no placebo
patient_data <- patient_data %>% 
  filter(idnr %in% as.numeric(colnames(counts_filtered))) %>% 
  filter(treatment != 1) %>% 
  drop_na()

#Add a column with the change in clinical parameters after 6 months
patient_data <- patient_data %>% 
  mutate(d_fev_pred = fevnapredv3 - fevnapredv1,
         d_ccq = ccqtotalv3 - ccqtotalv1,
         d_rvtlc = rvtlcpredv3 - rvtlcpredv1)

#Final selection of patients
counts_filtered <- counts_filtered %>%
  select(which(colnames(counts_filtered) %in% patient_data$idnr))

#Remove the "x" from the column names of norm_counts
if (exists("norm_counts")){
  norm_counts <- norm_counts %>% 
    rename_with(~gsub("X", "", .x))
}

