# code base for general purpose

* annotation_data - folder with files necessary for functional annotation of genes:
  * pathDB.rds - collection of Reactome, Kegg and MySigDB 50 Hallmark pathways
  * gomatrix.31.08.23.rds - GO-gene binary matrix: entrez IDs in rows, GO terms in columns
  * go2gName.rds - a list of GO terms, each element contains gene symdols of all unique genes AND genes of child terms
  * go2ensemble.rds - a list of GO terms, each element contains ensemble IDs of all unique genes AND genes of child terms
* func_analysis.R - a collection of functions for functional annotation of a gene/protein set
* usefulRfunc.r - a collection of useful R functions for data analysis and visualisation, collected over time
* justGO_ROBERT.R - GO enrichment test function that allows to adress redunduncy of GO terms by clustering and choosing the most representative among them, modified from https://github.com/robertsehlke/SETHRO

