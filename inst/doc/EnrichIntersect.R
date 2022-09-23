## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=TRUE---------------------------------------------------------------
library("EnrichIntersect")

## ----enrichment_map, results="hide", fig.width=7, fig.height=5----------------
x <- cancers_drug_groups$score
custom.set <- cancers_drug_groups$custom.set
set.seed(123)
enrich <- enrichment(x, custom.set, permute.n = 100)

## ----sankey_diagram, fig.width=9, fig.height=7--------------------------------
data(cancers_genes_drugs, package = "EnrichIntersect")
intersectSankey(cancers_genes_drugs, step.names = c("Cancers","Genes","Drugs"))

