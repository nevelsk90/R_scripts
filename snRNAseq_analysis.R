c# ############################################### #
# ####### analysis of Nphs2 mut. snRNAseq ####### #
# ############################################### #
library(ggplot2)
library(ggthemes)
library(plyr)
library(Matrix)
library(Seurat)
library(biomaRt)
library(SingleCellExperiment)

source("PROJECTS/myCode/Read10X_STAR.r")
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl", host="www.ensembl.org")
tx2gene <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),  mart = ensembl)



#### read in the data, create Seurat, do QC
  {
  matrix_dir = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_WT1hetdel/STAR/A006850158/GeneFull"
  ll <- list.dirs(  recursive = F , path=matrix_dir)
  annot_tab <- read.table(sep = "\t",header = T, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_WT1hetdel/STAR/A006850158/Sample_Names_A006850158.tab")
  # order annotation same as sequence files
  annot_tab <- annot_tab[ match( sub("_.*","",sub("^.*_1","1",ll)) , annot_tab$CCG.Sample.ID),]
  
  ### create seurat objects for all replicates
  # lapply( ll , function(x)
  #   {
    ### read exonic
    # Lane 1
    # y <- sub("_L001","_L002",x)
    
    filepath_premRNA <- paste0( ll, "/filtered/")
    names(filepath_premRNA) <- paste( annot_tab$CCG.Sample.ID , 
                                       paste(annot_tab$age , annot_tab$gtype, sep = "_") , sep = "__")
    

    # names( filepath_list ) <- c("Exon_L1","Exon_L2","Intr_L1","Intr_L2")
    snRNAseq_Wt1_premRNA <- Read10X_STAR( data.dir = filepath_premRNA , gene.column =2,
                                    zipped=F )
    
    
    # snRNAseq_Nphs2_Seurat <- CreateSeuratObject( counts =  snRNAseq_Nphs2_premRNA , project = "Nphs2_premRNA", min.cells = 100, min.features = 500 )
    sce <- CreateSeuratObject( counts = snRNAseq_Wt1_premRNA , project = "Wt1het.del_premRNA", min.cells = 100, min.features = 500 )
    # snRNAseq_Nphs2_Seurat <- CreateSeuratObject(counts = snRNAseq_Nphs2_intr , project = "Nphs2_intr", min.cells = 100, min.features = 500 )

  ### add idents and metadata
    sce$age <- ("4w")
    sce$gtype <- sub( ".*4_", "",sub( "wt_.*", "wt",
            colnames(sce)))
   
  ### QC 
  library(Seurat)
  # compute mt RNA
  sce$percent.mt <- PercentageFeatureSet(object = sce, pattern = "^mt-")

  # make QC plots
  # sce@meta.data$groups <- sub("Podocin.","",sce@meta.data$groups)
  cowplot::plot_grid( nrow = 1, rel_widths = c(1, -0, 1,-0,1),align = "hv",
                      VlnPlot(sce, features = c("nFeature_RNA"), 
                              pt.size = 0 , group.by = "orig.ident") +  guides(fill=FALSE) +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()) ,
                      NULL, VlnPlot(sce, features = c( "nCount_RNA" ), 
                                    pt.size = 0 , group.by = "orig.ident") + guides(fill=FALSE)+ 
                        theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()),
                      NULL,VlnPlot(sce, features = c( "percent.mt"), 
                                   pt.size = 0 , group.by = "orig.ident") + guides(fill=FALSE)+ 
                        theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank())
  ) 
  # filter cell with more than 1% of mt RNA
  sce4w <- subset( sce, cells = WhichCells(sce4w, expression = percent.mt < 1 ) )
  sce <- readRDS(  file= "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_WT1hetdel/Seurat/snRNAseq_Wt1hd_premRNA.Seurat.rda" )
  sce <- merge( sce4w, sce)
  saveRDS( sce, file= "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_WT1hetdel/Seurat/snRNAseq_Wt1hd_premRNA.4w12w25w.Seurat.rda" )
  
  }


#### do the analysis 
  {
  sce <- readRDS(  file= "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_WT1hetdel/Seurat/snRNAseq_Wt1hd_premRNA.Seurat.rda" )

  ### normialise the data
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  ### find variable features
  sce  <- FindVariableFeatures( sce , selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(sce), 10)   # Identify the 10 most highly variable genes
  
  ### scale the data
  sce <- ScaleData(sce)

  ### run PCA
  sce <- RunPCA( sce , features = VariableFeatures( object = sce ) )
  ElbowPlot( sce , ndims = 50 ) # choosing number of PCs to include 
  # DimPlot( sce, reduction = "pca", group.by = "gtype", order =  sample(colnames(sce)) )
  
  # cluster
  sce <- FindNeighbors(sce, dims = 1:20)
  sce <- FindClusters(sce, resolution = 0.1)
  
  # run UMAP
  sce <- RunUMAP(sce, dims = 1:20 )
  sce$age <- gsub( ".*__", "", as.factor( sce$age ))
  # saveRDS( sce , file= "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2_premRNA.Seurat.Norm.rda" )
  
  # plot of subsampled cells
  sce <- readRDS( file= "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2_premRNA.Seurat.Norm.rda" )
  sce_subs <- sce
  sce_subs$clust <- Idents(sce_subs) 
  Idents(sce_subs) <- sce@meta.data$groups
  sce_subs <- subset( sce_subs ,  downsample = 500  ) 
  Idents(sce_subs) <- sce_subs@meta.data$clust
  p1 <- DimPlot( sce_subs, reduction = "umap", order =  sample(colnames(sce_subs)) , label = T)
  p2 <- DimPlot( sce_subs, reduction = "umap", group.by = "groups", order =  sample(colnames(sce_subs)) )
  p3 <- DimPlot( sce_subs, reduction = "umap", group.by = "age", order =  sample(colnames(sce_subs)) )
  p4 <- DimPlot( sce_subs, reduction = "umap", group.by = "gtype", order =  sample(colnames(sce_subs)) )
  p5 <-   FeaturePlot( sce_subs, reduction = "umap", order = F,features = "nCount_RNA" , 
                 max.cutoff = 15000)
  p6 <-   FeaturePlot( sce_subs, reduction = "umap", order = F, features = "Wt1" , 
                       min.cutoff = 1)
  ### UMAPS to help with problematic clusters
  # Pax8, Dkk3, Igfbp6 (all parietal epithelial cell marker genes - cluster 13)
  # Ehd3, Ehd4 (may help to separate glomerular from extraglomerular endothelium and reveal more structure in cluster 0)
  # Aqp2 (marker gene for principal cells, will probably show up in cluster 4)
  # Slc8a1 (marker gene for connecting tubule, will probably show up in cluster 4, as well)
  FeaturePlot(sce_subs, cols=c("lightgrey", "red"),features = c("Pax8", "Dkk3", "Igfbp6","Ehd3","Ehd4","Aqp2","Slc8a1") )
  
  
  # podocyte marks
  VlnPlot(sce, features = c("Wt1", "Nphs1", "Magi2"), pt.size = 0 )
  
  
  # # PDS marks by groups
  # VlnPlot(sce, features = c("Wt1", "Mafb", "Foxd1"), 
  #         pt.size = 0 , group.by = "groups", idents =c("2","5"))
 VlnPlot(sce , features = c("Nphs1", "Nphs2","Wt1", "Mafb"), pt.size = 0 , 
          group.by = "groups", idents =c("2","5") , ncol = 4)
  
  
  # Podocyte markers
  # https://www.rngdsystems.com/resources/cell-markers/kidney-cells/podocytes/podocyte-markers
  cowplot::plot_grid(p1,p2, p3, p4,p5,p6, ncol = 3)
  # saveRDS(sce_subs, file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2_exon.Seurat.SubsTOviz.rda")
  sce_subs <- readRDS( file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2_premRNA.Seurat.SubsTOviz.rda")
  
  ### find markers
  sce_subsMarks <- FindAllMarkers( sce_subs )
  saveRDS( sce_subsMarks, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_Nphs2/Seurat/clustMarks_exon.rda")
  top10 <- sce_subsMarks %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  sce_subs_subs <- sce
  sce_subs_subs <- subset( sce_subs_subs ,  downsample = 200  ) 
  
  sce_subs_subs <- ScaleData(sce_subs_subs, features = rownames(sce_subs_subs) )
  DoHeatmap( sce_subs_subs, features = top10$gene) + NoLegend()
  
  ### compute doublets https://github.com/chris-mcginnis-ucsf/DoubletFinder
  {
    library(DoubletFinder) 
    library(ggplot2)
    ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    sweep.res.list_kidney <- paramSweep_v3(sce_subs, PCs = 1:30, sct = FALSE)
    sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
    bcmvn_kidney <- find.pK(sweep.stats_kidney) # pK=2.8
    
    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    annotations <- sce_subs@meta.data$clust
    homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
    nExp_poi <- round(0.05*nrow(sce_subs@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    # 
    sce_subs <- doubletFinder_v3(sce_subs, PCs = 1:30, pN = 0.25, pK = 0.28, nExp = nExp_poi.adj, 
                                 reuse.pANN = FALSE ,sct = FALSE )
    p1 <- DimPlot( sce_subs, reduction = "umap", order = F,
                   group.by  ="DF.classifications_0.25_0.28_341")
    p2 <- FeaturePlot( sce_subs, reduction = "umap", order = F,features = "nCount_RNA" , 
                       max.cutoff = 15000)
    cowplot::plot_grid(p1,p2, ncol = 2)
  }
  
  ### find doublets with scDblFinder https://bioconductor.org/books/release/OSCA/doublet-detection.html 
  {
    library(scDblFinder)
    library(SingleCellExperiment)
    # Doublet detection with clusters
    sce <- readRDS( file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2_premRNA.Seurat.Norm.rda")
    # sce_subs <- subset(sce, downsample=500)
    sce_scExp <- as.SingleCellExperiment( sce )
    dbl.out <- findDoubletClusters(sce_scExp, clusters=sce_scExp$seurat_clusters)
    dbl.out
    #  Doublet detection by simulation
    set.seed(100)
    
    # Setting up the parameters for consistency with denoisePCA();
    # this can be changed depending on your feature selection scheme.
    dbl.dens <- computeDoubletDensity(sce_scExp, d=ncol(reducedDim(sce_scExp)))
    summary(dbl.dens)
    
    sce_scExp$DoubletScore <- dbl.dens
    sce_scExp$DoubletScore_limit <- dbl.dens
    sce_scExp$DoubletScore_limit[sce_scExp$DoubletScore_limit > 5] <- 5
    
    scater::plotColData(sce_scExp, x="seurat_clusters", y="DoubletScore_limit", colour_by="seurat_clusters",
                        point_size =0.3)
 
    scater::plotUMAP(sce_scExp, colour_by="DoubletScore_limit", point_size =0.3  )
    scater::plotUMAP(sce_scExp, text_by="seurat_clusters" , colour_by="seurat_clusters", point_size =1   )
    
    # filter out doublets
    sce$DoubletScore <- dbl.dens
    saveRDS( sce , file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2_premRNA.Seurat.Norm.rda")
    
    sce_noDubl <- subset(sce, idents = c(0,1,2,3,4,5,7,8,10,12,13,14,15),
                         subset= DoubletScore <=0.5)
    scater::plotUMAP( as.SingleCellExperiment(sce_noDubl), text_by="seurat_clusters" , colour_by="seurat_clusters", point_size =0.5   )
    
    saveRDS( sce_noDubl , file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2_noDublets.12.04.21.rda")
    
    sce_noDub.Podo <- subset(sce_noDubl, idents=c(2,5))
    saveRDS( sce_noDub.Podo , file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2_Podo.NuDubs.rda")
    
     }
   
  ### visualise podocytes only
  {
    sce_noDub.Podo <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2_Podo.NuDubs.rda")
    ### find variable features
    sce_noDub.Podo  <- FindVariableFeatures( sce_noDub.Podo , selection.method = "vst", nfeatures = 2000)
    top10 <- head(VariableFeatures(sce_noDub.Podo), 10)   # Identify the 10 most highly variable genes
    
    ### scale the data
    sce_noDub.Podo <- ScaleData(sce_noDub.Podo)
    
    ### run PCA
    sce_noDub.Podo <- RunPCA(sce_noDub.Podo, features = VariableFeatures(object = sce_noDub.Podo))
    ElbowPlot(sce_noDub.Podo , ndims = 50 ) # choosing number of PCs to include 
    # DimPlot( sce, reduction = "pca", group.by = "gtype", order =  sample(colnames(sce)) )
    
    # cluster
    sce_noDub.Podo <- FindNeighbors(sce_noDub.Podo, dims = 1:10)
    sce_noDub.Podo <- FindClusters(sce_noDub.Podo, resolution = 0.05)
    sce_noDub.Podo$PDS  <- datPlot$score # add PDS score
    
    # run UMAP
    sce_noDub.Podo <- RunUMAP(sce_noDub.Podo, dims = 1:10 )
    saveRDS(readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2_Podo.NuDubs.rda"))
    
    ### make plots
    sce_noDub.Podo_SCE <- as.SingleCellExperiment(sce_noDub.Podo)
    set.seed(42)
    sce_noDub.Podo_subs <- subset( sce_noDub.Podo, cells=colnames(sce_subsPodo))
    sce_noDub.Podo_SCE <- as.SingleCellExperiment(sce_noDub.Podo_subs)
    # sce_noDub.Podo_SCE$PDS <- datPlot$score
    sce_noDub.Podo_SCE$PDS_lim <- datPlot$score
    sce_noDub.Podo_SCE$PDS_lim[sce_noDub.Podo_SCE$PDS_lim>-0.1] <- -0.1
    sce_noDub.Podo_SCE$PDS_lim[sce_noDub.Podo_SCE$PDS_lim<=-0.3] <- -0.3
    
    p1 <- scater::plotUMAP(sce_noDub.Podo_SCE, text_by="seurat_clusters" , colour_by="seurat_clusters", point_size = 2   )
    p2 <- scater::plotUMAP(sce_noDub.Podo_SCE, text_by="seurat_clusters" , colour_by="age", point_size = 2   )
    p3 <- scater::plotUMAP(sce_noDub.Podo_SCE, text_by="seurat_clusters" , colour_by="gtype", point_size = 2   )
    p4 <- scater::plotUMAP(sce_noDub.Podo_SCE, text_by="seurat_clusters" , colour_by="groups", point_size = 1   )
    p5 <- scater::plotUMAP(sce_noDub.Podo_SCE, text_by="seurat_clusters" , colour_by="PDS_lim", point_size = 2  )
    
    cowplot::plot_grid(p2, p3,p4, p5, ncol = 2)
    
    
    
    ### cluster markers
    sce_subsMarks <- FindAllMarkers( sce_noDub.Podo_subs )
    # saveRDS( sce_subsMarks, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_Nphs2/Seurat/clustMarks_exon.rda")
    top20 <- sce_subsMarks %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
    sce_subs_subs <- sce_noDub.Podo_subs
    sce_subs_subs <- subset( sce_subs_subs ,  downsample = 200  ) 
    
    sce_subs_subs <- ScaleData(sce_subs_subs, features = rownames(sce_subs_subs) )
    DoHeatmap( sce_subs_subs, features = top20$gene) + NoLegend()
    
    ### DE
    ages <- levels( as.factor(sce_noDub.Podo$age))
    Idents(sce_noDub.Podo) <- sce_noDub.Podo$groups
    DE.podo <- lapply( ages , function(xx) FindMarkers(sce_noDub.Podo ,
                                ident.2 = paste(xx,"Podocin.wt.wt",sep = "_") ,
                                ident.1 = paste(xx,"Podocin.R231Q.A286V",sep = "_")
                                   ))
    names(DE.podo) <- ages
    
    DE.podo_tab <- Reduce(rbind, lapply(seq(DE.podo), function(x) {
      DE.podo[[x]]$age <- ages[x]
      DE.podo[[x]]$gName <- rownames(DE.podo[[x]])
      return(DE.podo[[x]])
      
      ## venn diagram
      gplots::venn( lapply(DE.podo, rownames) )
      
    }))
    
  }
  
  ### trajectory for sick podocytes
  {
    ### plot cluster1
    scater::plotUMAP(sce_noDub.Podo_SCE, text_by="seurat_clusters" , colour_by="groups", point_size = 1   )
    sce_noDub.Podo_clust1 <- subset(sce_noDub.Podo, ident=1,subset= gtype=="Podocin.R231Q.A286V")
    sce_noDub.Podo_clust1SCE <- as.SingleCellExperiment(sce_noDub.Podo_clust1)
    pp1 <- scater::plotUMAP(sce_noDub.Podo_clust1SCE, colour_by="groups", point_size = 1   )
    pp3 <- scater::plotUMAP(sce_noDub.Podo_clust1SCE, colour_by="DoubletScore", point_size = 1   )
    pp4 <- scater::plotUMAP(sce_noDub.Podo_clust1SCE, colour_by="PDS", point_size = 1   )
    pp5 <- scater::plotUMAP(sce_noDub.Podo_clust1SCE, colour_by="nFeature_RNA", point_size = 1   )
    pp6 <- scater::plotUMAP(sce_noDub.Podo_clust1SCE, colour_by="Tmsb4x", point_size = 1   )
    pp7 <- scater::plotUMAP(sce_noDub.Podo_clust1SCE, colour_by="Spp1", point_size = 1   )
    pp8 <- scater::plotUMAP(sce_noDub.Podo_clust1SCE, colour_by="Epha6", point_size = 1   )
    
    pp9 <- scater::plotUMAP(sce_noDub.Podo_clust1SCE, colour_by="Rbfox1", point_size = 1   )
    pp10 <- scater::plotUMAP(sce_noDub.Podo_clust1SCE, colour_by="Igfbp5", point_size = 1   )
    pp11 <- scater::plotUMAP(sce_noDub.Podo_SCE, colour_by="Ephb1", point_size = 1   )
    pp12 <- scater::plotUMAP(sce_noDub.Podo_SCE, colour_by="Efnb1", point_size = 1   )
    pp13 <- scater::plotUMAP(sce_noDub.Podo_clust1SCE, colour_by="Ephb1", point_size = 2   )
    pp14 <- scater::plotUMAP(sce_noDub.Podo_clust1SCE, colour_by="Efnb1", point_size = 2   )
    
    sce_noDub.Podo_clust1SCE$Efnb1_limit <- sce_noDub.Podo_clust1SCE@assays@data$logcounts["Efnb1",]
    sce_noDub.Podo_clust1SCE$Efnb1_limit[sce_noDub.Podo_clust1SCE$Efnb1_limit >2.5] <- 2.5
    sce_noDub.Podo_clust1SCE$Ephb1_limit <- sce_noDub.Podo_clust1SCE@assays@data$logcounts["Ephb1",]
    sce_noDub.Podo_clust1SCE$Ephb1_limit[sce_noDub.Podo_clust1SCE$Ephb1_limit >2.5] <- 2.5
    
    #### https://bioconductor.org/books/release/OSCA/trajectory-analysis.html
    library(scater)
    pseudo.all <- TSCAN::quickPseudotime(sce_noDub.Podo_clust1SCE,
                                         clusters = sce_noDub.Podo_clust1SCE$age, 
                                         use.dimred="UMAP",  with.mnn=F)
    mnn.pseudo <- rowMeans(pseudo.all$ordering, na.rm=TRUE)
    pp2 <- scater::plotUMAP(sce_noDub.Podo_clust1SCE, colour_by=I(mnn.pseudo), text_by="age", text_colour="red") +
      geom_line(data=pseudo.all$connected$UMAP, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))
    
    # Changes along a trajectory
    pseudo <- TSCAN::testPseudotime(sce_noDub.Podo_clust1SCE , pseudotime=mnn.pseudo )
    pseudo <- pseudo[order(pseudo$p.value),]
    
    cowplot::plot_grid(pp1 , pp2)
    cowplot::plot_grid( pp4, pp5, pp6 , pp7,
                       pp8,pp9,  ncol = 3)
    cowplot::plot_grid(pp11 , pp12)
    cowplot::plot_grid(pp13 , pp14)
    
    }
    
  }

#### DE and trajectory analysis for all clusters ####
  {
    ### load subsampled data, no doublets removal
    load( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Sniny_snRNAseq_Nphs2mut/snRNAseqNphs2mut_shiny.Rdata")
    
    kidney_Nphs2mut_Seurat$Cell_Type <- factor( kidney_Nphs2mut_Seurat$seurat_clusters, levels = levels(Idents(kidney_Nphs2mut_Seurat)), 
                                                labels = c("endothelium","mesangium","podocytes","proximal tubules",
                                                           "(distal) connecting tubules, CDPC","podocytes","doublets: endothelium + mesangium",
                                                           "loop of Henle (ascending)", "immune cells","doublets: mesangium + podocytes",
                                                           "loop of Henle (descending)","doublets: endothelium + podocytes", 
                                                           "juxtaglomerular apparatus", "S3 proximal tubule","collecting duct intercalated cells",
                                                           "cycling cells", "doublets") )
    Idents(kidney_Nphs2mut_Seurat) <- kidney_Nphs2mut_Seurat$Cell_Type
    
    ###
  kidney_Nphs2mut_ctrVSmut.DE <- lapply( levels(kidney_Nphs2mut_Seurat$Cell_Type),
                                         function(ii) {
                                           
                                           w4_DE <- tryCatch( FindMarkers( kidney_Nphs2mut_Seurat, subset.ident = ii , group.by = 'groups',
                                                                           ident.1 = "4-5_Podocin.wt.wt", ident.2 = "4-5_Podocin.R231Q.A286V" ,
                                                                           min.cells.group	= 50 ), error = function(e){ return(NULL) } )
                                           if( !is.null(w4_DE) ) { 
                                             w4_DE$Gene.name <- rownames(w4_DE)
                                             w4_DE$age <- rep("4", nrow(w4_DE))
                                           } 
                                           
                                           w6_DE <- tryCatch( FindMarkers( kidney_Nphs2mut_Seurat, subset.ident = ii , group.by = 'groups',
                                                                           ident.1 ="6_Podocin.wt.wt", ident.2 = "6_Podocin.R231Q.A286V" ,
                                                                           min.cells.group	= 50 ), error = function(e){ return(NULL) } )
                                           if( !is.null(w6_DE) ) { 
                                             w6_DE$Gene.name <- rownames(w6_DE)
                                             w6_DE$age <- rep("6", nrow(w6_DE))
                                           }
                                           
                                           w8_DE <- tryCatch( FindMarkers( kidney_Nphs2mut_Seurat, subset.ident = ii , group.by = 'groups', 
                                                                           ident.1 ="8_Podocin.wt.wt", ident.2 = "8_Podocin.R231Q.A286V",
                                                                           min.cells.group	= 50 ), error = function(e){ return(NULL) } )
                                           if( !is.null(w8_DE) ) { 
                                             w8_DE$Gene.name <- rownames(w8_DE)
                                             w8_DE$age <- rep("8", nrow(w8_DE))}
                                           
                                           w12_DE <- tryCatch( FindMarkers( kidney_Nphs2mut_Seurat, subset.ident = ii , group.by = 'groups',
                                                                            ident.1 = "12_Podocin.wt.wt", ident.2 = "12_Podocin.R231Q.A286V" ,
                                                                            min.cells.group	= 50 ), error = function(e){ return(NULL) } )
                                           if( !is.null(w12_DE) ) { 
                                             w12_DE$Gene.name <- rownames(w12_DE)
                                             w12_DE$age <- rep("12", nrow(w12_DE))
                                           }
                                           
                                           
                                           DEtab <- rbind(w4_DE , w6_DE, w8_DE, w12_DE)
                                           return(DEtab)
                                           
                                         } )
  
  library(scater)
  kidney_Nphs2mut_mutProg.DE <- lapply( levels(kidney_Nphs2mut_Seurat$Cell_Type), 
                                        function(ii){
                                          sce_clust <- subset( kidney_Nphs2mut_Seurat , ident=ii, subset= gtype=="Podocin.R231Q.A286V")
                                          sce_clust_SCE <- as.SingleCellExperiment(sce_clust)
                                          
                                          pseudo.all <- TSCAN::quickPseudotime( sce_clust_SCE,
                                                                                clusters = sce_clust_SCE$age, 
                                                                                use.dimred="UMAP",  with.mnn=F)
                                          
                                          mnn.pseudo <- rowMeans(pseudo.all$ordering, na.rm=TRUE)
                                          
                                          ffile <- paste("/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/snRNAseq_podo/Seurat/premRNA/DE/",
                                                         "Nphs2_FSGStraj_clust_",ii, ".pdf",sep = "")
                                          cairo_pdf( ffile, height=5, width=10, fallback_resolution=1200 )
                                          
                                          pp1 <- scater::plotUMAP(sce_clust_SCE, colour_by="groups", point_size = 1   )
                                          pp2 <- scater::plotUMAP(sce_clust_SCE, colour_by=I(mnn.pseudo), text_by="age", text_colour="red") +
                                            geom_line(data=pseudo.all$connected$UMAP, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))
                                          
                                          
                                          
                                          cowplot::plot_grid(pp1 , pp2)
                                          dev.off()
                                          
                                          # Changes along a trajectory
                                          pseudo <- TSCAN::testPseudotime( sce_clust_SCE , pseudotime=mnn.pseudo )
                                          pseudo <- pseudo[order(pseudo$p.value),]	
                                          
                                          return(pseudo)
                                          
                                        })
  
}

### calculate PDS 
  {
  sce <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2_noDublets.12.04.21.rda")
  # ### extract cluster  2,5 as clear podocyte clusters
  sce_subsPodo <- subset( sce , ident=c("2","5"), downsample= 3000)
  # sce_subsPodo <- subset( sce_subs , ident=c("6"))
  # sce_noDub.Podo <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2_")
# 
#     # make annotation
#     datt_des <- paste( sce$age ,sce$gtype, sep = "_")
#     datt_des <- relevel( as.factor( datt_des ), ref = "4-5_Podocin.wt.wt" )

  # load final list of marker genes
  meanRANK_marks_all <- readRDS( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/FSGS_markers0.75_36stud.rds" )
  # get gene sets
  geneSets <- meanRANK_marks_all[[2]]
  
  # calculate for each cluster seperately
  PDS.50_tab <- Reduce( rbind, lapply( levels(Idents( sce ) ),
                                       function( clust){
                                         datt <- subset(sce , idents = clust)
                                         
                                         # saveRDS( snRNAseq_podo@assays$RNA@data , file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_WT1hetdel/Seurat/snRNAseq_Wt1hd_Podo.07.04.21.rda")
                                         # rownames( snRNAseq_podo ) <- toupper( rownames( snRNAseq_podo ) )
                                         # snRNAseq_podo <- cbind.data.frame( Gene.symbol =rownames( snRNAseq_podo ), snRNAseq_podo )
                                         
                                         
                                         
                                         # 1. Build gene-expression rankings for each cell  
                                         exprMat <- as.matrix( datt@assays$RNA@data )
                                         rownames( exprMat ) <- toupper( rownames( datt ) )
                                         cells_rankings <-  AUCell::AUCell_buildRankings( exprMat , nCores=1, plotStats=F)
                                         
                                         
                                         # 2. Calculate enrichment for the gene signatures (AUC)
                                         cells_AUC <- AUCell::AUCell_calcAUC(geneSets, cells_rankings )
                                         cells_AUC <- AUCell::getAUC( cells_AUC)
                                         
                                         # data for plotting
                                         design_podo <- as.factor( paste(datt$age, datt$gtype,sep = "_"))
                                         
                                         datPlot <- cbind.data.frame( score= ( cells_AUC[4,]-cells_AUC[5,] ),
                                                                      genotype_age = design_podo )
                                         
                                         return(datPlot)
                                         gc()
                                       }) )
  
  # PDS.50_tab
  saveRDS( PDS.50_tab, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/snRNAseq_Nphs2/PDS/snRNAseq_NoDubs_PDStop50_Nphs2.rds")
    # PDS.50_tab <- readRDS()
  
  samples <- sce$groups
  # modify group names to exclude sample IDs
  sce$groups <- relevel( as.factor( paste( sce$age , sce$gtype,sep = "_") ) , 
                         ref = "4-5_Podocin.wt.wt" )
  # remove duplicated cluster column
  sce@meta.data <- sce@meta.data[, -which(names(sce@meta.data) %in% "RNA_snn_res.0.1")]
  # add PDS as metadata
  sce$PDS.50 <- PDS.50_tab$score[ match( colnames(sce), rownames(PDS.50_tab))]
  saveRDS( sce , file= "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_WT1hetdel/Seurat/snRNAseq_Wt1hd_noDublets.PDS50.19.04.21.rda" )
  
  ### gene set with podomarkers only
  {
    sce_subsMarks <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_Nphs2/Seurat/clustMarks_premRNA.rda")
    sce_subsMarks_podo <- unique( sce_subsMarks$gene[which(sce_subsMarks$cluster%in%c(2,5) & 
                                                             sce_subsMarks$p_val_adj<0.01 & 
                                                             abs(sce_subsMarks$avg_logFC)>1)] )
    
    geneSets_podoMark <-  GSEABase::GeneSetCollection( c(
      GSEABase::GeneSet( geneSets[[4]]@geneIds[geneSets[[4]]@geneIds %in% toupper(sce_subsMarks_podo)], setName="top50 UP markers") ,
      GSEABase::GeneSet( geneSets[[5]]@geneIds[geneSets[[5]]@geneIds %in% toupper(sce_subsMarks_podo)] , setName="top50 DOWN markers"),
      GSEABase::GeneSet( geneSets[[8]]@geneIds[geneSets[[8]]@geneIds %in% toupper(sce_subsMarks_podo)] , setName="top100 UP markers") , 
      GSEABase::GeneSet( geneSets[[9]]@geneIds[geneSets[[9]]@geneIds %in% toupper(sce_subsMarks_podo) ] , setName="top100 DOWN markers")
    ))
    
    # 2. Calculate enrichment for the gene signatures (AUC)
    cells_AUC <- AUCell::AUCell_calcAUC(geneSets_podoMark, cells_rankings )
    cells_AUC <- AUCell::getAUC( cells_AUC)
    
    # data for plotting
    datPlot <- cbind.data.frame( score= ( cells_AUC[3,]-cells_AUC[4,] ),
                                 genotype_age = sce_subsPodo_des )
    
    # plot 
    ggplot2::ggplot( datPlot , aes( x=score, color=genotype_age)) + 
      scale_color_colorblind() +  geom_density(size=1.5) + theme_bw() + 
      labs(x = "podocyte damage score" , title = "snRNA-seq Nphs2 clusters 2, 5") +
      geom_vline(data=ddply( datPlot, "genotype_age", summarise, grp.mean=mean(score)), aes(xintercept=grp.mean, color=genotype_age), linetype="dashed")
    
    }
  
 }

