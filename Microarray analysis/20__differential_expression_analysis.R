#Select relevant columns from patient data
clin_data <- patient_data %>% select(patientID, starts_with("d_"))

for (i in 2:ncol(clin_data)){
  print(colnames(clin_data[i]))
  clin_param <- clin_data %>% pull(i)
  
  #Design matrix
  mm <- model.matrix(~ clin_param)
  row.names(mm) <- clin_data$patientID
  
  #voom transformation
  d <- voom(expression_data_cleaned)
  
  #Estimate the correlation between duplicate spots
  corfit <- duplicateCorrelation(d, mm)
  
  #Fit linear model
  fit <- lmFit(d, mm, correlation = corfit$consensus.correlation)
  
  efit <- eBayes(fit)
  
  #Get results
  tt <-topTable(efit,
                adjust = "BH",
                number = nrow(efit))
  
  file_name <- paste0("DEG_", colnames(clin_data[i]), "_v3_baseline", ".csv")
  write.table(x = tt, 
              file = here("Microarray analysis", "output", file_name), 
              sep = ",")
}
