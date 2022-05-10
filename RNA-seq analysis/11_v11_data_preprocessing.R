#Filter samples available at baseline
counts_filtered_v11 <- counts %>% 
  select(ends_with("_1")) %>% 
  rename_with(~gsub("_1|X", "", .x))

#Select required columns
patient_data_v11 <- patient_data_raw %>% 
  select(idnr, treatment, geslacht, lftv1, rookjaar, rooknuv1, packyearv1, 
         ccqtotalv1, ccqtotalv11, rvtlcpredv1, rvtlcpredv11, fevnapredv1, fevnapredv11) 

#Filter for patients that have RNA-seq available at baseline, no NA in any column and have received no placebo
patient_data_v11 <- patient_data_v11 %>% 
  filter(idnr %in% as.numeric(colnames(counts_filtered_v11))) %>% 
  filter(treatment != 1) %>% 
  filter(treatment != 4) %>%
  drop_na()

#Add a column with the change in clinical parameters after 6 months
patient_data_v11 <- patient_data_v11 %>% 
  mutate(d_fev_pred = fevnapredv11 - fevnapredv1,
         d_ccq = ccqtotalv11 - ccqtotalv1,
         d_rvtlc = rvtlcpredv11 - rvtlcpredv1)

#Final selection of patients
counts_filtered_v11 <- counts_filtered_v11 %>%
  select(which(colnames(counts_filtered_v11) %in% patient_data_v11$idnr))

#Get CPM normalized counts
norm_counts_v11 <- counts_filtered_v11 %>% 
  cpm(log = T) %>% 
  as.data.frame()


# #Get normalized counts
# ##Create a DGEList object
# d0 <- DGEList(counts_filtered_v11)
# 
# ##Filter lowly expressed genes
# keep <- filterByExpr(d0)
# d0 <- d0[keep, , keep.lib.sizes = F]
# 
# ##Normalize counts
# d1 <- calcNormFactors(d0)
# 
# ##Get TMM normalized counts
# norm_counts_v11 <- as.data.frame(cpm(d1, log = T))
