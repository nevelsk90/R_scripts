# code base for general purpose

* annotation_data - folder with files necessary for functional annotation of genes:
  * pathDB.rds - collection of Reactome, Kegg and MySigDB 50 Hallmark pathways
  * gomatrix.31.08.23.rds - GO-gene binary matrix: entrez IDs in rows, GO terms in columns
  * go2gName.rds - a list of GO terms, each element contains gene symdols of all unique genes AND genes of child terms
  * go2ensemble.rds - a list of GO terms, each element contains ensemble IDs of all unique genes AND genes of child terms
* func_analysis.R - a collection of functions for functional annotation of a gene/protein set
* usefulRfunc.r - a collection of useful R functions for data analysis and visualisation, collected over time
* justGO_ROBERT.R - GO enrichment test function that allows to adress redunduncy of GO terms by clustering and choosing the most representative among them, modified from https://github.com/robertsehlke/SETHRO

## Single-cell / spatial helpers
* Read10X_STAR.r - load STARsolo / 10X count matrices (modified Seurat `Read10X`)
* Create_10X_H5_starSolo.r - write STARsolo output as a 10X-style HDF5 file (modified from Seurat, GPL-3)
* stacked_vlnPlot_seurat.R - stacked violin plots for Seurat objects
* split_viol_func.r - split violin geom for ggplot2 (from a StackOverflow answer, link in the file)
* scVar.R - identify biological variability in scRNA-seq data (author: Huiwen Zheng)
* crop_seurat_v1.R, SeuratCropper.R - crop tissue areas / FOVs of spatial Seurat objects
* TFA.estimate_functions.r - transcription factor activity from expression and a TF-target prior (generalised inverse)

## Pathways, networks, statistics
* SPIA_analysis.r - signalling pathway impact analysis (SPIA via graphite)
* Helper_chordPlot.R, chord_plot.r - circular chord plots of genes and pathways
* layout_modular.r - igraph layout that groups nodes by community (modified from StackOverflow)
* extendGenes.coord_script.r - extend gene coordinates (GRanges)
* rbo_tim.R - rank-biased overlap between ranked lists (author: Fabian Schmich, gespeR)
* stabpath_diffPenalty.R - stability paths with different penalties per predictor
* induction_repression_functions.R - induction / repression time points in time series
* example_module.R - minimal Shiny module example

## Shell / pipelines
* star_bulk_rnaseq.sh - STAR index + mapping for bulk RNA-seq
* QuantSalmon_script.sh - Salmon quantification from a docker image
* Genrich_script.sh - ChIP/ATAC peak calling with Genrich
* get_SRR_data.sh - download SRA runs listed in SRR_Acc_List.txt
* nfcore_pipelines/ - SLURM scripts for nf-core/fetchngs and nf-core/rnaseq
* Pipelines/ - BWA mapping + QC for ChIP-seq (SLURM and plain shell)
* uni_script.sh - collection of command-line snippets (SRA, QC, MACS2, bowtie, ...); facility password via `$SEQ_PASSWORD`
