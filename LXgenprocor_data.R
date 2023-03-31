
#setwd("D:/Desktop/R包开发/LXDESeq")
library(openxlsx)

gene_data_example <- read.xlsx("gene_data.xlsx")
protein_data_example <- read.xlsx("protein_data.xlsx")


usethis::use_data(gene_data_example,overwrite = T)
usethis::use_data(protein_data_example,overwrite = T)

rm(list=ls())

data(gene_data_example)
data(protein_data_example)

