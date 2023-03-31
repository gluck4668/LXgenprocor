\name{LXgenprocor}
\alias{LXgenprocor}
\title{LXgenprocor}
\usage{
LXgenprocor(gene_data,protein_data,FC,species)
}
\description{
LXgenprocor can be used to analyze the correlation between the differently expressed genes (DEGs) and the differently expressed proteins (DEPs).
}
\examples{

if(!requireNamespace("devtools"))
  install.packages("devtools")

library(devtools)

install_github("gluck4668/LXgenprocor")

library(LXgenprocor)

??LXgenprocor

#------------------------------------
data(gene_data_example)
data(protein_data_example)

#------------------------------------

gene_data <- "gene_data.xlsx"

protein_data <- "protein_data.xlsx"

FC <- NULL # 1 or NULL

species <- " mouse" # Here, it should be "human", "mouse", or "rat"


LXgenprocor(gene_data,protein_data,FC,species)


}
