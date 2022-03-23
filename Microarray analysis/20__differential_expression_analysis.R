#Design matrix
mm <- model.matrix(~ patient_data$d_fev_pred)
row.names(mm) <- patient_data$patientID

#voom transformation
d <- voom(expression_data_cleaned, mm)

#Estimate the correlation between duplicate spots
corfit <- duplicateCorrelation(d, mm)




fit <- lmFit(d, mm, correlation = corfit$consensus.correlation)

efit <- eBayes(fit)

tt <-topTable(efit)
tt

https://rdrr.io/bioc/limma/man/dupcor.html
https://here.r-lib.org/