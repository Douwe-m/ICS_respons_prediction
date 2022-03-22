#-------------Perform DGE analysis-------------
#read data files
source("RNA-seq analysis/00_load_data.R", echo = T)

#Pre-process data
source("RNA-seq analysis/10_data_preprocessing.R", echo = T)

#Get DEGs
source("RNA-seq analysis/20_gene_expression_analysis.R", echo = T)


#-------------Visualize RNA-seq results-------------
#read data files
source("RNA-seq analysis/00_load_data.R", echo = T)

#Pre-process data
source("RNA-seq analysis/10_data_preprocessing.R", echo = T)

#Visualize the results
source("RNA-seq analysis/30_data_visualisation.R", echo = T)


#-------------GSVA analysis-------------
#read data files
source("RNA-seq analysis/00_load_data.R", echo = T)

#Pre-process data
source("RNA-seq analysis/10_data_preprocessing.R", echo = T)

#GSVA analysis
source("RNA-seq analysis/40_GSVA_analysis.R", echo = T)


#-------------DGE analysis over 30 months-------------
#read data files
source("RNA-seq analysis/00_load_data.R", echo = T)

#Pre-process data and DGE analysis
source("RNA-seq analysis/21_v11_vs_baseline.R")

for (i in res_v11){
  print(i[1:5,])
}

#-------------DGE analysis with outliers removed-------------
#Patients 2109 and 2107 were removed since they are possible outliers 
#read data files
source("RNA-seq analysis/00_load_data.R", echo = T)

#Pre-process data
source("RNA-seq analysis/10_data_preprocessing.R", echo = T)

#Perform analysis again, but with outliers (2109 and 2107) removed
source("RNA-seq analysis/22_gene_expression_analysis_patients_removed.R", echo = T)

for (i in top_genes) {
  print(i[1:5,])
}

#-------------DGEs with FDR <= 0.10 compared to clin. para. 6 months-------------
#read data files
source("RNA-seq analysis/00_load_data.R", echo = T)

#Pre-process data
source("RNA-seq analysis/10_data_preprocessing.R", echo = T)

DEG_d_ccq %>% 
  filter(FDR <= 0.10) %>% 
  pull(ensembl_gene_id)

DEG_d_fev %>% 
  filter(FDR <= 0.10) %>% 
  pull(ensembl_gene_id)

DEG_d_rvtlc %>% 
  filter(FDR <= 0.10) %>% 
  pull(ensembl_gene_id)

#GSVA analysis with change in clinical parameters after 6 months
source("RNA-seq analysis/40_GSVA_analysis.R", echo = T)

#-------------DGEs with FDR <= 0.10 compared to clin. para. 30 months-------------
#read data files
source("RNA-seq analysis/00_load_data.R", echo = T)

#Pre-process data
source("RNA-seq analysis/11_v11_data_preprocessing.R", echo = T)

#GSVA analysis with change in clinical parameters after 30 months
source("RNA-seq analysis/40_GSVA_analysis.R", echo = T)
