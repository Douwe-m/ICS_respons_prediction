#read count data
counts <- read.delim("D:/AFO/Data/GLUCOLD/GLUCOLD - mRNA Gene Expression (RNA Seq)/GLUCOLDhumangeneexpression.txt", 
                     sep = ";", header = T, row.names = 1)

#Gene symbols
gene_symbols <- read.delim("D:/AFO/Data/GLUCOLD/GLUCOLD - mRNA Gene Expression (RNA Seq)/genes_info_hg37.csv", 
                           sep = ",", header = T)

#load patient data
patient_data <- haven::read_sav("D:/AFO/Data/GLUCOLD/GLUCOLD - Patients/2015-11-24 GLUCOLD database (totaal, 114pt).sav")