#Design matrix
mm <- model.matrix(~ patient_data$d_fev_pred)
row.names(mm) <- patient_data$patientID

fit <- lmFit(expression_data_cleaned, mm)

efit <- eBayes(fit)

tt <-topTable(fit)
tt
