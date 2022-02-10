# ############################################### #
# ####### analysis of Nphs2 mut. snRNAseq ####### #
# ############################################### #
library(ggplot2)
library(ggthemes)
library(plyr)
library(dplyr)
library(Matrix)
library(Seurat)
library(biomaRt)
library(SingleCellExperiment)
library(matrixStats)

source("PROJECTS/myCode/Read10X_STAR.r")
source( "/home/tim_nevelsk/PROJECTS/myCode/usefulRfunc.r")

ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl", host="www.ensembl.org")
tx2gene <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),  mart = ensembl)



#### read new 6-8w samples, create Seurat, do QC, merge with other samples
{
  
  matrix_dir = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/STAR/STARpremRNA/A006850158"
  ll <- list.dirs(  recursive = F , path=matrix_dir)
  annot_tab <- read.table(sep = "\t",header = T, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/STAR/STARpremRNA/A006850158/Sample_Names.tab")
  # order annotation same as sequence files
  annot_tab <- annot_tab[ match( sub("_.*","",sub("^.*_1","1",ll)) , annot_tab$CCG.Sample.ID),]
  
 
  filepath_premRNA <- paste0( ll, "/filtered/")
  names(filepath_premRNA) <- paste( annot_tab$CCG.Sample.ID , 
                                    paste(annot_tab$age , annot_tab$gtype, sep = "_") , sep = "__")
  
  
  # names( filepath_list ) <- c("Exon_L1","Exon_L2","Intr_L1","Intr_L2")
  snRNAseq_Nphs2_premRNA <- Read10X_STAR( data.dir = filepath_premRNA , gene.column =2,
                                        zipped=F )
  
  
  # snRNAseq_Nphs2_Seurat <- CreateSeuratObject( counts =  snRNAseq_Nphs2_premRNA , project = "Nphs2_premRNA", min.cells = 100, min.features = 500 )
  sce <- CreateSeuratObject( counts = snRNAseq_Nphs2_premRNA , project = "Nphs2mut_premRNA", min.cells = 100, min.features = 500 )
  # snRNAseq_Nphs2_Seurat <- CreateSeuratObject(counts = snRNAseq_Nphs2_intr , project = "Nphs2_intr", min.cells = 100, min.features = 500 )
  
  ### add idents and metadata
  sce$age <-  sub( "_mut.*", "",sub( ".*__", "", colnames(sce)))
  sce$gtype <- sub( ".*mut", "mut",sub( "mut_.*", "mut",
                                    colnames(sce)))
  

  ### combine with the previous samples
  sce <- readRDS(  file= "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2_premRNA.Seurat.rda" )
  sce <- merge( sce6w, sce)
  sce$gtype <- as.character( factor( sub( "_.*" , "", sce$gtype) , 
                                     levels = c( "mut" ,"Podocin.R231Q.A286V" ,"Podocin.wt.wt" ),
                                     labels = c("Nphs2mut","Nphs2mut","wtype") )) 
  sce$groups <- paste(sce$gtype, sce$age , sep = "_")

  ### QC 
  # compute mt RNA
  sce$percent.mt <- PercentageFeatureSet(object = sce, pattern = "^mt-")
  # filter cell with more than 1% of mt RNA
  sce <- subset( sce, cells = WhichCells( sce , expression = percent.mt < 1 ) )
  # make QC plots
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
  
  saveRDS( sce, file= "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2mut.allSamples_premRNA.Seurat.rda" )
  
  
}


#### do the analysis 
{
  sce <- readRDS(  file= "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2mut.allSamples_premRNA.Seurat.rda" )
  
  ### normialise the data
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  ### find variable features
  sce <- FindVariableFeatures( sce , selection.method = "vst", nfeatures = 2000)

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
  sce$groups <-paste( sce$gtype, sce$age , sep = "_")
  
  saveRDS( sce, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2mut.allSamples_premRNA.Seurat.Norm.rda")
  
  # plot of subsampled cells
  set.seed(42)
  sce_subs <- subset( sce ,  downsample = 500  ) 
  
  p1 <- DimPlot( sce_subs, reduction = "umap", shuffle = T , label = T)
  p2 <- DimPlot( sce_subs, reduction = "umap", group.by = "groups", shuffle = T  )
  p3 <- DimPlot( sce_subs, reduction = "umap", group.by = "age", shuffle = T  )
  p4 <- DimPlot( sce_subs, reduction = "umap", group.by = "gtype", shuffle = T  )
  p5 <- FeaturePlot( sce_subs, reduction = "umap", order = F,features = "nCount_RNA" , 
                     max.cutoff = 15000)
  p6 <-   FeaturePlot( sce_subs, reduction = "umap", order = F, features = "Wt1" , 
                       min.cutoff = 0)
  # combine plots
  pdf( height = 12 , width = 18, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/snRNAseq_Nphs2/Seurat/FinalKFO/snRNAseq_Nphs2mut.allSamples_UMAPs.pdf")
   cowplot::plot_grid( p1,p2, p3, p4,p5,p6, ncol = 3)
  dev.off()
  
  ### Violin plot for podocyte marks
  VlnPlot(sce, features = c("Wt1", "Nphs1", "Magi2","Mafb","Foxd1","Plekhh2"), pt.size = 0 )  
  
  ### do PCA on pseudobulk
  sce_subs <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2mut.allSamples_premRNA.Seurat.Viz.rda")
  XX <- as.data.frame(t(sce_subs@assays$RNA@counts))
  XX$groups <- as.factor(sce_subs$orig.ident)
  XX <-   aggregate( .~ groups , data= XX, sum )
  datPlot <- t(XX[,-1])
  colnames(datPlot) <- XX[,1]
  # normalise by library size 
  datPlot <- t(t(datPlot) / colSums(datPlot))*1000000
  
    # plot contribution of genes to PCs 
  factoextra::fviz_cos2( prcomp(t(log(datPlot+1)),scale=F), 
                         choice = "var", axes = 1:2, top=20)

  
  PCA_heatm_plot( count_table = datPlot , 
                  sampleNames= levels( as.factor(paste(sce_subs$orig.ident, sce_subs$gtype, sce_subs$age, sep = "_")))  , 
                  groupsExp = c( "wt", "mut", "wt", "wt", "mut", "mut", "mut",
                                 "wt","mut","wt","wt", "mut","mut","mut") , 
                  logtrans= T , title="PCA on pseudobulk Nphs2 snRNAseq",
                  topNprct=0.1 )
  
  ### find markers
  sce_subsMarks <- FindAllMarkers( sce_subs )
  saveRDS( sce_subsMarks, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2_premRNA.allSamples.clustMarks.rda")
  top10 <- sce_subsMarks %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  # downsample the data  further for a better heatmap
  sce_subs_subs <- subset( sce_subs ,  downsample = 200  ) 
  # scale the data for all genes to avoid missing genes on the heatmap
  sce_subs_subs <- ScaleData(sce_subs_subs, features = rownames(sce_subs_subs) )
  DoHeatmap( sce_subs_subs, features = top10$gene) + NoLegend()
  
  ### extract podocytes
  sce_Podo <- subset( sce , idents=c( 2 , 5 ))
  # saveRDS( sce_Podo, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2mut.allSamples_premRNA.Seurat.Podo.rda")
  
  ### compute doublets https://github.com/chris-mcginnis-ucsf/DoubletFinder
  {
    library(DoubletFinder) 
    library(ggplot2)
    sce_subs <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2mut.allSamples_premRNA.Seurat.Viz.rda")
    ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    sweep.res.list_kidney <- paramSweep_v3(sce_subs, PCs = 1:20, sct = FALSE)
    sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
    par(mar=c(5,5,2,2))
    bcmvn_kidney <- find.pK(sweep.stats_kidney) # pK=0.26
    
    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    annotations <- sce_subs@meta.data$clust
    homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
    nExp_poi <- round(0.1*nrow(sce_subs@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    # 
    sce_subs <- doubletFinder_v3(sce_subs, PCs = 1:20, pN = 0.25, 
                                 pK = 0.26, nExp = nExp_poi.adj, 
                                 reuse.pANN = FALSE ,sct = FALSE )
    p1 <- DimPlot( sce_subs, reduction = "umap", order = F,
                   group.by  ="DF.classifications_0.25_0.26_800")
    p2 <- FeaturePlot( sce_subs, reduction = "umap", order = F,features = "nCount_RNA" , 
                       max.cutoff = 15000)
    cowplot::plot_grid(p1,p2, ncol = 2)
  }
  
  ### find doublets with scDblFinder https://bioconductor.org/books/release/OSCA/doublet-detection.html 
  {
    library(scDblFinder)
    library(SingleCellExperiment)
    # Doublet detection with clusters
    # sce <- readRDS( file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2_premRNA.Seurat.Norm.rda")
    # sce_subs <- subset(sce, downsample=500)
    sce_scExp <- as.SingleCellExperiment( sce_subs )
    dbl.out <- findDoubletClusters(sce_scExp, 
                clusters=sce_scExp$seurat_clusters , 
                get.all.pairs = T )
    dbl.out
    
    #  Doublet detection by simulation
    set.seed(42)
    
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
}

### analyse podocytes 
{
  sce_Podo <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2mut.allSamples_premRNA.Seurat.Podo.rda" )
  set.seed(42)
  sce_subsPodo <- subset( sce_Podo , downsample=2000)
  # Idents(sce_subsPodo) <- sce_subsPodo$orig.ident
  # sce_subsPodo <- subset( sce_subsPodo , downsample=500)
  # Idents(sce_subsPodo) <- sce_subsPodo$seurat_clusters
  
  ### PCA with podocytes
  {
    XX <- as.data.frame( t( as.matrix( sce_subsPodo@assays$RNA@counts) ) )
    XX$groups <- as.factor( sce_subsPodo$orig.ident)
    XX <-   aggregate( .~ groups , data= XX, sum )
    datPlot <- t(XX[,-1])
    colnames(datPlot) <- XX[,1]
    # normalise by library size 
    datPlot <- t(t(datPlot) / colSums(datPlot))*1000000
    
    # # plot contribution of genes to PCs 
    # factoextra::fviz_cos2( prcomp(t(log(datPlot+1)),scale=F), 
    #                        choice = "var", axes = 1:2, top=20)
    # 
    
    PCA_heatm_plot( count_table = datPlot , 
                    sampleNames= levels( as.factor(paste(sce_subsPodo$orig.ident, sce_subsPodo$gtype, sce_subsPodo$age, sep = "_")))  , 
                    groupsExp = c( "wt", "mut", "wt", "wt", "mut", "mut", "mut",
                                   "wt","mut","wt","wt", "mut","mut","mut") , 
                    logtrans= T , title="PCA on pseudobulk Nphs2 snRNAseq Podocytes",
                    topNprct=0.1 )
  }

  ### visualise podocytes
  {
    ### find variable features
    sce_Podo <- FindVariableFeatures( sce_Podo , selection.method = "vst", nfeatures = 2000)
    
    ### scale the data
    sce_Podo <- ScaleData(sce_Podo)
    
    ### run PCA
    sce_Podo <- RunPCA( sce_Podo , features = VariableFeatures( object = sce_Podo ) )
    ElbowPlot( sce_Podo , ndims = 50 ) # choosing number of PCs to include 
    
    # cluster
    sce_Podo <- FindNeighbors( sce_Podo, dims = 1:10 )
    sce_Podo <- FindClusters( sce_Podo, resolution = 0.05 )
    
    # run UMAP
    sce_Podo <- RunUMAP( sce_Podo , dims = 1:10 )
    
    ### make plots
    sce_subsPodo <- subset( sce_Podo , downsample=2000)
    sce_Podo_SCE <- as.SingleCellExperiment(sce_subsPodo)
    p1 <- scater::plotUMAP(sce_Podo_SCE, text_by="seurat_clusters" , colour_by="seurat_clusters", point_size = 2   )
    p2 <- scater::plotUMAP(sce_Podo_SCE, text_by="seurat_clusters" , colour_by="age", point_size = 2   )
    p3 <- scater::plotUMAP(sce_Podo_SCE, text_by="seurat_clusters" , colour_by="gtype", point_size = 2   )
    p4 <- scater::plotUMAP(sce_Podo_SCE, text_by="seurat_clusters" , colour_by="groups", point_size = 1   )
    p5 <- scater::plotUMAP(sce_Podo_SCE, text_by="seurat_clusters" , colour_by="PDS.50.Podo", point_size = 2  )
    p6 <- scater::plotUMAP(sce_Podo_SCE, text_by="seurat_clusters" , colour_by="Epha6", point_size = 2  )
    
    cowplot::plot_grid(p1, p2, p3,p4,p5,p6, ncol = 3)
  }
  
  ### do DE for podocytes
  {

    # check clusters that look like doublets
    # 
    DE.podo <- list( 
      {
        X <- FindMarkers(sce_subs , ident.2 =  5, ident.1 = 2 )
        X$test <- "clust5_VS_clust2" 
        X$gNames <- rownames(X)
        X
      } ,
      {
        X2 <- FindMarkers(sce_subs ,ident.2 = 9,ident.1 = 2 )
        X2$test <- "clust9_VS_clust2" 
        X2$gNames <- rownames(X2)
        X2
      }, 
      {
        X3 <-FindMarkers(sce_subs ,ident.2 = 12 ,ident.1 = 2) 
        X3$test <- "clust12_VS_clust2" 
        X3$gNames <- rownames(X3)
        X3
      }
    ) 
    names(DE.podo) <- c("clust5_VS_clust2",  
                        "clust9_VS_clust2",
                        "clust12_VS_clust2")
    # save the table
    DE.podo_tab <- Reduce( rbind, DE.podo)
    DE.podo_tab$status <- ifelse( DE.podo_tab$gNames %in% Reduce( intersect, lapply(DE.podo, function(X) X$gNames)), 
                                  "DE_inAll", ifelse( DE.podo_tab$gNames %in% intersect(DE.podo[[2]]$gNames, DE.podo[[3]]$gNames) ,
                                                      "DE_double",  "") )
    DE.podo_tab$sig <- ifelse( DE.podo_tab$p_val_adj < 0.1, "YES","NO")
    write.table(DE.podo_tab , sep = "\t", quote = F, row.names = F, 
                file="/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/snRNAseq_Nphs2/Seurat/FinalKFO/snRNAseq_Nphs2_DE.double.tsv")
    
    ## venn diagram
    gplots::venn( lapply(DE.podo, rownames) )
    
    
  }
  
  ### trajectory for sick  podocytes
  {
    ### plot cluster1
    # sce_PodoClust2 <- subset( sce_subsPodo, idents=c(  2 , 5) )
    sce_PodoClust5 <- subset( sce_Podo, idents=c( 5), subset= gtype=="Nphs2mut" )
    
     sce_Podo_SCE <- as.SingleCellExperiment( sce_PodoClust5 )
     sce_Podo_SCE <- scater::arrange( data = sce_Podo_SCE , sample( colnames(sce_Podo_SCE)))
     
    scater::plotUMAP( sce_Podo_SCE , text_by="seurat_clusters" , colour_by="groups", point_size = 2   )
    
    #### https://bioconductor.org/books/release/OSCA/trajectory-analysis.html
    library(scater)
    pseudo.all <- TSCAN::quickPseudotime( sce_Podo_SCE ,
                                          clusters = sce_Podo_SCE$groups , 
                                          use.dimred = "UMAP" , with.mnn=F )
    
    mnn.pseudo <- rowMeans( pseudo.all$ordering , na.rm=TRUE )
    pp1 <-  scater::plotUMAP(sce_Podo_SCE, colour_by="groups", point_size = 1   ) + 
      xlim(NA,-8) + ylim(-2,3)
    
    pp2 <- scater::plotUMAP(sce_Podo_SCE, colour_by=I(mnn.pseudo), text_by="groups", 
                            text_colour="red", point_size = 1 ) + xlim(NA,-8)+ ylim(-2,3)+
      geom_line(data=pseudo.all$connected$UMAP, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))
    pp3 <- scater::plotUMAP(sce_Podo_SCE, colour_by="PDS.50.Podo", text_by="groups", 
                            text_colour="red", point_size = 1 ) + xlim(NA,-8)+ ylim(-2,3)+
      geom_line(data=pseudo.all$connected$UMAP, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))
    pp4 <- scater::plotUMAP(sce_Podo_SCE, colour_by="Epha6", text_by="groups", 
                            text_colour="red", point_size = 1 ) + xlim(NA,-8)+ ylim(-2,3)+
      geom_line(data=pseudo.all$connected$UMAP, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))
    # pdf( height= )
    cowplot::plot_grid(pp1 , pp2 , pp3,pp4)
    
    # # Changes along a trajectory
    pseudo <- TSCAN::testPseudotime(sce_Podo_SCE , pseudotime=mnn.pseudo )
    pseudo <- pseudo[order(pseudo$p.value),]
    
    
    
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

#### calculate PDS ####
{
  sce <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2_noDublets.12.04.21.rda")

  ### calculate for podocytes
  {
    # load final list of marker genes
    meanRANK_marks_all <- readRDS( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/FSGS_markers0.75_36stud.rds" )
    # get gene sets
    geneSets <- meanRANK_marks_all[[2]]
    
    ### calculate PDSfor podocytes
    # 1. Build gene-expression rankings for each cell  
    exprMat <- as.matrix( sce_Podo@assays$RNA@data )
    rownames( exprMat ) <- toupper( rownames( sce_Podo ) )
    cells_rankings <-  AUCell::AUCell_buildRankings( exprMat , nCores=1, plotStats=F)
    
    
    # 2. Calculate enrichment for the gene signatures (AUC)
    cells_AUC <- AUCell::AUCell_calcAUC(geneSets, cells_rankings ,
                                        aucMaxRank = ceiling(0.2 * nrow(cells_rankings)))
    cells_AUC <- AUCell::getAUC( cells_AUC)
    
    # data for plotting
    sce_Podo$PDS.50.Podo <- ( cells_AUC[4,]-cells_AUC[5,] )
    sce_Podo$PDS.20.Podo <- ( cells_AUC[2,]-cells_AUC[3,] )
    sce_Podo$PDS.200.Podo <- ( cells_AUC[6,]-cells_AUC[7,] )
    
    ### plot distributions 
    # sce_subsPodo <- sce_Podo
    # Idents(sce_subsPodo) <- sce_subsPodo$groups
    # sce_subsPodo <- subset( sce_subsPodo , downsample=1000)
    XX <- sce_Podo@meta.data
    # XX <- XX[XX$seurat_clusters %in% c(12),]
    ggplot2::ggplot( XX , aes( x=PDS.200.Podo, color= groups)) + 
      scale_color_colorblind() +  geom_density(size=1.5) + theme_bw() + 
      labs(x = "podocyte damage score" , title = "snRNA-seq Nphs2 podocytes") +
      theme( text = element_text( size = 20)) +
      geom_vline( data=ddply( XX , "groups", summarise, 
                              grp.mean=mean(PDS.200.Podo)), aes(xintercept=grp.mean, 
                                                                color=groups), linetype="dashed" )
    ###plot using vaseplots
    # test differnce between means
    my_comparisons <- compare_means(formula = PDS.50.Podo ~ gtype ,  
                                    data = XX  ) 
    
    # make a vase plot https://github.com/heike/ggboxplots/blob/master/R/stat-hdr.r
    library(ggboxplots)
    library( ggpubr )
    
    na.rm=T
    ggplot2::ggplot( XX , aes( y=PDS.50.Podo, 
      x = reorder( orig.ident , PDS.50.Podo, FUN = median) , color=gtype) ) +
      geom_vase(fill = "white", width = 1,lwd=1.5) + 
      scale_color_colorblind() + theme_bw() +
      # coord_cartesian(ylim = quantile( sceHomoPodo@meta.data$PDS.50, c(0.02, 0.98) ))+
      theme( text =  element_text(size=20) ,
             axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      labs( title="Nphs2 mut. snRNAseq data", x ="samples ID, ordered by median" ,
            subtitle= paste( "Wilcoxon test: ",my_comparisons$group1[1],
                             "VS",my_comparisons$group2[1],"p =",my_comparisons$p.format[1])) 
    
  }

  ### calculate for each cluster seperately
  {
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
    
  }
  
}

