#Design matrix
mm <- model.matrix(~ patient_data$d_ccq)
row.names(mm) <- patient_data$patientID

d <- voom(expression_data_cleaned, mm)

fit <- lmFit(d, mm)

efit <- eBayes(fit)

tt <-topTable(efit)
tt

https://rdrr.io/bioc/limma/man/dupcor.html
https://here.r-lib.org/