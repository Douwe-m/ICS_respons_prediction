#-------------Perform DGE analysis-------------
#read data files
source("RNA-seq analysis/00_load_data.R", echo = T)

#Pre-process data
source("RNA-seq analysis/10_data_preprocessing.R", echo = T)

#Get DEGs
source("RNA-seq analysis/20_gene_expression_analysis.R", echo = T)


#-------------Identified DEGs with change in clin. para. after 6 months -------------
#read data files
source("RNA-seq analysis/00_load_data.R", echo = T)

#Pre-process data
source("RNA-seq analysis/10_data_preprocessing.R", echo = T)
source("RNA-seq analysis/11_v11_data_preprocessing.R", echo = T)


#GSVA analysis
source("RNA-seq analysis/40_GSVA_analysis.R", echo = T)


#-------------Identified DEGs with change in clin. para. after 30 months -------------
#read data files
source("RNA-seq analysis/00_load_data.R", echo = T)

#Pre-process data
source("RNA-seq analysis/11_v11_data_preprocessing.R", echo = T)

#GSVA analysis
source("RNA-seq analysis/40_GSVA_analysis.R", echo = T)

