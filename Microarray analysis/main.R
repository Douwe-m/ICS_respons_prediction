#-------------Main analysis-------------
#read data files
source("Microarray analysis/00__load_data.R", echo = T)

#Pre-process data
source("Microarray analysis/10__data_preprocessing.R", echo = T)

#Get DEGs
source("Microarray analysis/20__differential_expression_analysis.R", echo = T)

#-------------Insepect DEGs-------------
#read data files
source("Microarray analysis/00__load_data.R", echo = T)

#Pre-process data
source("Microarray analysis/10__data_preprocessing.R", echo = T)

#Inspect DEGs
source(here::here("Microarray analysis", "30__gsva_analysis.R"), echo = T)

#-------------GSVA analysis with change in clin param over 30 months-------------
#read data files
source(here::here("Microarray analysis", "00__load_data.R"), echo = T)

#Pre-process data
source(here::here("Microarray analysis", "11__v11_data_preprocessing.R"), echo = T)

#GSVA analysis
source(here::here("Microarray analysis", "30__gsva_analysis.R"), echo = T)
