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

#Visualize the results
source("RNA-seq analysis/30_data_visualisation.R", echo = T)

#GSVA analysis
source("RNA-seq analysis/40_GSVA_analysis.R", echo = T)

#Prepare data for MLR
source(here::here("RNA-seq analysis", "50_MLR_data_prep.R"), echo = T)

#Feature selection and MLR
source(here::here("RNA-seq analysis", "52_MLR_analysis.R"), echo = T)



#-------------Identified DEGs with change in clin. para. after 30 months -------------
#read data files
source("RNA-seq analysis/00_load_data.R", echo = T)

#Pre-process data
source("RNA-seq analysis/11_v11_data_preprocessing.R", echo = T)

#Visualize the results
source("RNA-seq analysis/30_data_visualisation.R", echo = T)

#GSVA analysis
source("RNA-seq analysis/40_GSVA_analysis.R", echo = T)



#-------------Other analyses-------------
#-------------DGE analysis over 30 months-------------
#read data files
source("RNA-seq analysis/00_load_data.R", echo = T)

#Pre-process data and DGE analysis
source("RNA-seq analysis/21_v11_vs_baseline.R")

#-------------DGE analysis with outliers removed-------------
#Patients 2109 and 2107 were removed since they are possible outliers 
#read data files
source("RNA-seq analysis/00_load_data.R", echo = T)

#Pre-process data
source("RNA-seq analysis/10_data_preprocessing.R", echo = T)

#Perform analysis again, but with outliers (2109 and 2107) removed
source("RNA-seq analysis/22_gene_expression_analysis_patients_removed.R", echo = T)
