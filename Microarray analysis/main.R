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
source(here::here("Microarray analysis", "30__inspect_deg.R"), echo = T)
