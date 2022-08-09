#### analysis of bulk RNAseq from KFO ####
options(connectionObserver = NULL)

library(ggplot2)
library(biomaRt)
library(viridis)
library(reshape2)
library( ggthemes)


### load DE results with shrunk lfc (for meanLFC and gsea analysis)
ll <- list("/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_Wt1/DEexonIntron_4w_glvl.lfcShrunk.rda",
           "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_Wt1/DEexonIntron_12w_glvl.lfcShrunk.rda",
           "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_NPHS2/Nphs2_IntExon_glvl_4wDE.lfcShrunk.rda",
           "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_NPHS2/Nphs2_IntExon_glvl_8wDE.lfcShrunk.rda", 
           "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_adriamycin/adriamycin.bulkRNAseq_glvl_DE.lfcShrunk.rda")

DElist <- lapply( ll, readRDS )
names(DElist) <- c( "Wt1het.del_4w","Wt1het.del_12w","Nphs2.mut_4w","Nphs2.mut_8w","adriamycin_4w" )

####  compare global results #### 
# heatmap of LFC
ggname <- na.omit(Reduce(union, lapply(DElist, function(X) rownames(X))))
LFC_tab <- Reduce(cbind, lapply(DElist, function(X) X$lfc_shrunk[match(ggname, rownames(X))]))
colnames(LFC_tab) <- names(DElist)
rownames(LFC_tab) <- ggname
LFC_tab <- LFC_tab[rowSums(is.na(LFC_tab))<4,]
# LFC_tab[is.na(LFC_tab)] <- 0
corrplot::corrplot(cor(LFC_tab, use="pairwise.complete.obs"))

### LFC distribution
datToplot <- LFC_tab[,1:4]
datToplot <- melt(datToplot )
datToplot[is.na(datToplot)] <- 0
datToplot$experiment <- sub("_.*", "", datToplot$Var2)
datToplot$stage <-  sub(".*_", "", datToplot$Var2)
datToplot$stage <- ifelse(datToplot$stage=="4w", "early ", "late") 

ggplot( data = datToplot, aes(x=(value), color=experiment, linetype=stage))+ 
  geom_density( lwd=1.5, bw=0.05, na.rm = T) + 
  theme_bw()+ scale_colour_grey()+ xlab("log2FC")+
  theme( text = element_text(size = 20)) + labs(color='Datasets') + 
  coord_cartesian(xlim = c(-0.5,0.5),ylim = c(0,6) ) 
  

#### intersect of DE genes
library(UpSetR)
ll <- lapply(1:4, function(ii){ rownames(DElist[[ii]])[
  DElist[[ii]]$padj<0.05 & !is.na(DElist[[ii]]$padj)]})
names(ll)<- names(DElist)[1:4]
upset( fromList(ll ), point.size = 3, line.size = 1.5, 
       number.angles = 30, text.scale = 2)

#### perform GO and pathway annotation ####
### use Robert function
{
  options(connectionObserver = NULL)
  source("/home/tim_nevelsk/PROJECTS/myCode/justGO_ROBERT.R")  
  library("org.Mm.eg.db")
  
  # use Robert's function to create sparse binary matrix indicating membership of genes (columns) in GO terms (rows)
  gomatrix=sf.createGoMatrix()
  # load biomart
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
  tx2gene_GO <- biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id"), mart = mart)
  
  # run function on list of DE results
  GOrobert_res <- lapply(seq(DElist), function(ii){
    print(ii)
    
    datt <- DElist[[ii]]
    # define background gene set
    universe= unique(tx2gene_GO$entrezgene_id[match(  rownames(datt), 
                                                   tx2gene_GO$ensembl_gene_id)]) 
    universe <- as.character( universe[!is.na(universe)] )
    
    
    # prepare gene set, convert ensembleIDs to entrezIDs
    geneset=rownames(datt)[ datt$padj<0.05 & !is.na(datt$padj) ]
    geneset <- unique(tx2gene$entrezgene_id[match( geneset,  tx2gene_GO$ensembl_gene_id)])
    geneset <- geneset[!is.na(geneset)]
    geneset <- as.character(geneset)
    print(length(intersect(geneset,colnames(gomatrix))))
    
    # apply Robert's function that given a sparse matrix of GO terms (columns = genes, rows = GO terms)
    # a geneset of interest and a background set of genes (universe)
    # carry out clustering with members diverging by at most cut_max genes, and do enrichment testing.
    # Note, multiplicity adjustment is performed for the representative terms only.
    RobertGO=sf.clusterGoByGeneset( gomatrix, geneset, universe, 
                                    min.genes=5, cut_max = 5 )
    
    return(RobertGO)
  })
  saveRDS(GOrobert_res, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_adriamycin/KFO.bulkRNAseq_GOrobert.rda")
  
  # make plots
  lapply( seq(GOrobert_res), function(ii){
    print(ii)
    datt <- GOrobert_res[[ii]]
    datt <- subset(datt$results, Enrichment>0 & Is.primary==TRUE & Primary.Fisher.adj<0.1)
    datt <- datt[ order(datt$Primary.Fisher.adj),][1:30, ]
    
    margin_x=c(0,12,2,0)
    pp<- sf.plotResults(  datt )
    print(pp)
  } )
  
}

### use topGO
{
  # prepare GO list
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "mmusculus_gene_ensembl",host = 'ensembl.org')
  All_GO=getBM(attributes=c("ensembl_gene_id","go_id"), mart=mart)
  GO_list=unstack(All_GO ) 
  GO_list <- GO_list[names(GO_list)!=""]
  
  topGO_res <- lapply( seq(DElist), function(ii){
    print(ii)
    
    datt <- DElist[[ii]]
    # define the universe
    universe <- as.character(  rownames(datt) )
    # prepare gene set, convert ensembleIDs to entrezIDs
    geneset=rownames(datt)[datt$padj<0.05 & !is.na(datt$padj) ]
    geneset <- as.character(geneset)
    print(length(intersect(geneset,unlist(GO_list))))
    
    GoBytopGO_weight01( GO_list=GO_list, geneset = geneset ,universe = universe , min.genes=3 )
    
  } )
  saveRDS(topGO_res , file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_adriamycin/KFO.bulkRNAseq_topGO.rda")

  
  }

### GSEA
{
  library(fgsea)
  
  # Reactome
  {
    # read the file
    reactPath <- fgsea::gmtPathways( "/media/tim_nevelsk/WD_tim/ANNOTATIONS/Pathways/Reactome/ReactomePathways.gmt")
    # convert human to mouse IDS
    
    fun_homoTOmouse <- function(gns){
      egs <- mapIds(org.Hs.eg.db, gns, "ENTREZID","SYMBOL")
      mapped <- select(Orthology.eg.db, egs, "Mus.musculus","Homo.sapiens")
      mapped$MUS <- mapIds(org.Mm.eg.db, as.character(mapped$Mus.musculus), "SYMBOL", "ENTREZID")
      mapped$MUS.ens <- mapIds(org.Mm.eg.db, as.character(mapped$Mus.musculus), "ENSEMBL", "ENTREZID")
    return(mapped)
      }
    
    options(connectionObserver = NULL)
    library( "org.Hs.eg.db" )
    library("org.Mm.eg.db")
    library(Orthology.eg.db)
    
    humanGen <- unique (Reduce( c , reactPath) )
    mouseGen <- fun_homoTOmouse(humanGen)
    mouseGen$HOMO <- rownames(mouseGen)
    
    reactPath_gName <- lapply( reactPath , function(x){
      X <- as.character( na.omit( unique(mouseGen$MUS.ens[ match( x, mouseGen$HOMO)])))
    } )
    
    # remove untranslated names
    reactPath_gName <- lapply(reactPath_gName, function(X) X <- X[X!="NULL"])
    
    # # delete paths with less than 3 genes
    # reactPath_gName <- reactPath_gName[ unlist( lapply( reactPath_gName, function(x) length(x)>2  ))]
    
    }
  
  # KEGG
  {
    ## download directly from KEGG: 341 paths
    library(EnrichmentBrowser)
    options(connectionObserver = NULL)
    
    library(KEGGgraph)
    keggPath_gName <- EnrichmentBrowser::getGenesets( "mmu", db ="kegg",
                                                      gene.id.type = "ENSEMBL" )
    
    names(keggPath_gName) <- names(keggPath_gName)
    
     
  }
  
  # use shrunk lfc
  fgsea_res <- lapply( seq(DElist), function(ii){
    print(ii)
    
    datt <- DElist[[ii]]
    datt <- setNames( datt[["lfc_shrunk"]] , rownames(datt))
    datt <- datt[order(-datt)]                              
    fgseaREACT <- fgseaMultilevel(pathways = reactPath_gName, 
                      stats    = datt,
                      minSize  = 5  ,maxSize = 500 )
    fgseaKEGG <- fgseaMultilevel(pathways = keggPath_gName, 
                        stats    = datt,
                        minSize  = 5 ,maxSize = 500 )
    ll <- list( fgseaKEGG=fgseaKEGG , fgseaREACT=fgseaREACT)
    return(ll)
  })
  
  saveRDS( lapply(fgsea_res, "[[" , 1) , file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/Analysis/bulk/KFO.bulkRNAseq_fgsea.KEGG.rda")
  saveRDS( lapply(fgsea_res, "[[" , 2) , file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/Analysis/bulk/KFO.bulkRNAseq_fgsea.REACT.rda")
  
  }



#### compare GO and pathway annotations ####
GOrobert_res <- readRDS( file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_adriamycin/KFO.bulkRNAseq_GOrobert.cut50.rda")
topGO_res <- readRDS( file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_adriamycin/KFO.bulkRNAseq_topGO.rda")
KEGG_res <- readRDS( file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_adriamycin/KFO.bulkRNAseq_fgsea.KEGG.rda")
REACT_res <-  readRDS( file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_adriamycin/KFO.bulkRNAseq_fgsea.REACT.rda")

### create a df with average LFC for all  GO terms
GO_meanLFC <- Reduce( rbind, lapply( seq(unique(All_GO$go_id)) , function(ii){
  # print(ii)
  datt <- unique(All_GO$go_id)[ii]
  ggenes <- All_GO$ensembl_gene_id[All_GO$go_id==datt ]
  if( !is.matrix(LFC_tab[rownames(LFC_tab)%in%ggenes,]) ) {
    # filter out GOs with less than 3 DE genes detected in bulkRNAseq
    c(0,0,0,0) } else if( nrow(LFC_tab[rownames(LFC_tab)%in%ggenes,]) < 3 ){
      c(0,0,0,0) }  else colMeans(LFC_tab[rownames(LFC_tab)%in%ggenes,], na.rm = T)
}))
rownames(GO_meanLFC) <- unique(All_GO$go_id)
GO_meanLFC <- GO_meanLFC[ rowSums(GO_meanLFC)!=0 & !is.na(rowSums(GO_meanLFC)) ,]
GO_meanLFC.4 <- GO_meanLFC[,1:4]
# remove GO terms with litle variation
hist(log(rowVars(GO_meanLFC.4)))
GO_meanLFC_sel <- GO_meanLFC.4[log(rowVars(GO_meanLFC.4)) > (mean(log(rowVars(GO_meanLFC.4))) +
                                                           2*sd(log(rowVars(GO_meanLFC.4)))), ]
  
### PCA analysis of GO functions
{
 
  ### do PCA
  library(factoextra)
  library(viridis)
  # match IDs to term
  library(GO.db)
  GO <- as.list(GOTERM)
  rownames( GO_meanLFC_sel ) <- sapply( GO[rownames(GO_meanLFC_sel)], 
                                        function(X) wrap_text( X@Term))
  # run PCA  
  res.pca <- prcomp( t(GO_meanLFC_sel[,1:4]), scale = F)
  
  # plot results
  fviz_eig(res.pca)
  
  # Contributions of variables to PC1
  a<-fviz_contrib(res.pca, choice = "var", axes = 1,top=10)
  # Contributions of variables to PC2
  b<-fviz_contrib(res.pca, choice = "var", axes = 2,top=10)
  # 
  gridExtra::grid.arrange(grobs= lapply(list(a,b), "+", theme(plot.margin=margin(10,10,2,2))), ncol=2, 
                          top='Contribution of the variables to the first two PCs')
  
  # GO terms that separate datasets based on LFC.PCA
  var <- get_pca_var(res.pca)
  GOdiffLFC <- Reduce( union , list( rownames(var$cos2[order(-var$cos2[,1]),][1:10,]),
                                     rownames(var$cos2[order(-var$cos2[,2]),][1:10,])))
  
  fviz_pca_biplot(res.pca, repel = TRUE, geom.ind = c("point", "text"),
                  geom.var = c("point", "text") ,  axes = c(1, 2),
                  select.var = list( name = GOdiffLFC ), 
                  labelsize = 6, 
                  col.var = "#2E9FDF", # Variables color
                  col.ind = "#696969"  # Individuals color
  )   + theme( text = element_text(size = 22 ),
               axis.text.x = element_text(size = 14))  
}


#### buble plot of GO and pathways ####
###  GO terms
{
  # get a set of common func
  GOrobert_common  <- Reduce( intersect, lapply( 1:4, function(ii, thrsh=0.1,
                                                               datt=GOrobert_res){
    print(ii)
    datt <- datt[[ii]]
    datt <- subset(datt$results,  Is.primary==TRUE & Primary.Fisher.adj<thrsh)
    datt$Term
  }))
  # # topGO
  # topGO_common <-  Reduce( intersect, lapply( 1:5, function(ii){
  #   print(ii)
  #   datt <- topGO_res[[ii]]
  #   datt <- subset(datt[,1:7], weight01<0.1)
  #   datt$Term
  # }))
  
  # get extra information from robertGO result table
  XX<- GOrobert_res[[1]]$results
  XX <- XX[ XX$Term%in%GOrobert_common , ]
  datTOplot <- GO_meanLFC[, 1:4]
  # datTOplot <- scale(datTOplot)
  datTOplot <- as.data.frame( datTOplot[rownames(datTOplot) %in%XX$GO.ID,] )
  datTOplot$size <- XX$Annotated[ match( rownames(datTOplot),XX$GO.ID)]
  datTOplot$padj <- XX$Fisher[ match( rownames(datTOplot),XX$GO.ID)]
  datTOplot$Ontology <- XX$Term[ match( rownames(datTOplot),XX$GO.ID)]
  datTOplot$dbase <- XX$Ontology[ match( rownames(datTOplot),XX$GO.ID)]
  # melt for ggplot
  datTOplot.melt1 <- reshape2::melt(datTOplot, id=c("size","padj","Ontology","dbase"))
  # plot
  ggplot(datTOplot.melt1, aes(x=value, y=Ontology,fill=-log2(padj),size=size)) + 
    scale_fill_distiller( direction = 1)+
    xlab("mean log2FC")+ scale_size(range = c(0, 15))+
    geom_vline(xintercept=0, linetype=2,lwd=2)+ # add vertical line for lfc=0
    geom_point(alpha=0.7,pch=21,stroke = 2) + facet_grid(cols = vars(variable), rows = vars(dbase),
                                       scales = "free_y", space = "free")+
    theme_bw() + theme( text = element_text(size = 22 ),
                        axis.text.x = element_text(size = 14)) 
  
}

# common pathways
{
  KEGG_res <- readRDS( file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/Analysis/bulk/KFO.bulkRNAseq_fgsea.KEGG.rda")[1:4]
  REACT_res <-  readRDS( file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/Analysis/bulk/KFO.bulkRNAseq_fgsea.REACT.rda")[1:4]

  # get top 10 pathways
  ppath <- Reduce( intersect, lapply(1:4, function(ii) KEGG_res[[ii]]$pathway ))
  KEGG_rank  <- rowMeans( Reduce( cbind, lapply( 1:4, function(ii){
    print(ii)
    datt <- KEGG_res[[ii]]
    datt <- datt[ match(ppath , datt$pathway),]
    setNames(  rank(datt$pval ), ppath)
  })))
  KEGG_rank <- KEGG_rank[order(KEGG_rank)]
  
  ppath <- Reduce( intersect, lapply(1:4, function(ii) REACT_res[[ii]]$pathway ))
  REACT_rank  <- rowMeans( Reduce( cbind, lapply( 1:4, function(ii){
    print(ii)
    datt <- REACT_res[[ii]]
    datt <- datt[ match(ppath , datt$pathway),]
    setNames(  rank(datt$pval ), ppath)
  })))
  REACT_rank <- REACT_rank[order(REACT_rank)]
  

  XX<- rbind( KEGG_res[[1]], REACT_res[[1]] )
  # get extra information from robertGO result table
  REACT_common <- names(REACT_rank)[1:20]
  KEGG_common <- NULL
  datTOplot <- rbind.data.frame( Reduce( cbind , lapply(KEGG_res, function(X) X$NES[match( KEGG_common, X$pathway )])) , 
                      Reduce( cbind , lapply(REACT_res, function(X) X$NES[match( REACT_common, X$pathway )])))
  rownames(datTOplot) <- c(KEGG_common, REACT_common)
  colnames(datTOplot) <- colnames(LFC_tab)[1:4]
  datTOplot$dbase <- c( rep("KEGG", length(KEGG_common)), 
                            rep("REACTOME",length(REACT_common)))
  datTOplot$size <- XX$size[ match( rownames(datTOplot),XX$pathway)]
  datTOplot$padj <- XX$pval[ match( rownames(datTOplot),XX$pathway)]
  datTOplot$Ontology <- XX$pathway[ match( rownames(datTOplot),XX$pathway)]
  # melt for ggplot
  datTOplot.melt2 <- reshape2::melt(datTOplot, id=c("size","padj","Ontology","dbase"))
  # datTOplot.melt <- rbind( datTOplot.melt1 , datTOplot.melt2)
  # plot
  ggplot(datTOplot.melt2, aes(x=value, y=Ontology ,fill=-log2(padj),size=size)) + 
    scale_fill_distiller( direction = 1)+xlab("NES of gsea")+scale_size(range = c(0, 15))+
    geom_vline(xintercept=0, linetype=2)+ # add vertical line for lfc=0
    geom_point(alpha=0.7,pch=21,stroke = 2) + 
     facet_grid(cols = vars(variable), rows = vars(dbase), scales = "free_y", space = "free")+
    theme_bw() + theme( text = element_text(size = 22 ),
                          axis.text.x = element_text(size = 14)) 
  

}

#### 2D comparison of annotation ####
### wt1
# get GO terms significant in at least one stage
X <- union(GOrobert_res[[1]]$results$GO.ID[which(GOrobert_res[[1]]$results$Fisher<0.01)],
           GOrobert_res[[2]]$results$GO.ID[which(GOrobert_res[[2]]$results$Fisher<0.01)] )
# get mean LFC for the terms of interest
GOrobert_LFC <- as.data.frame( GO_meanLFC[rownames(GO_meanLFC) %in% X,1:2] )
# add term size
GOrobert_LFC$size <- sapply( 1:nrow(GOrobert_LFC), function(ii){
  length( All_GO$ensembl_gene_id[All_GO$go_id == rownames(GOrobert_LFC)[ii]])
} )
# shorten the names
GOrobert_LFC$Primary <- sapply( GO[rownames(GOrobert_LFC)], 
                                function(X) wrap_text( X@Term))
# Create prediction interval data frame with upper and lower lines corresponding to sequence covering minimum and maximum of x values in original dataset
newx <- lm( GOrobert_LFC$Wt1het.del_12w ~ GOrobert_LFC$Wt1het.del_4w)
pred_interval <- predict(newx, interval="prediction", level = 0.95)
pred_interval <- as.data.frame(pred_interval)
GOrobert_LFC <- cbind(GOrobert_LFC , pred_interval)
# 
GOrobert_LFC$diffGO <- ifelse( GOrobert_LFC$Wt1het.del_12w < GOrobert_LFC$lwr |  
                                 GOrobert_LFC$Wt1het.del_12w > GOrobert_LFC$upr ,
                               as.character(GOrobert_LFC$Primary), "")

# plot
  library( ggrepel)

ggplot( GOrobert_LFC , aes(x= Wt1het.del_4w ,y= Wt1het.del_12w) ) +  
   geom_point(aes(size = size), colour="red",shape=21, stroke = 2) + 
  geom_line( data=GOrobert_LFC, aes(x = fit, y = upr), colour = "darkgrey")+
  geom_line( data=GOrobert_LFC, aes(x = fit, y = lwr), colour = "darkgrey")+
  scale_size(range = c(1, 15)) +  # adjust range of bubbles sizes
  geom_text_repel( aes(label = diffGO ), max.overlaps=100) + #ensure that labels are not jammed together
  xlab("mean log2FC week 4") + ylab("mean log2FC week 12") + 
  ggtitle("activity changes in GO terms") + theme_bw()+
  theme( text = element_text(size = 24))



#### comparison of PPI (STRING) ####
#### integration with TRN (ATACseq) ####
ATACseq_tgenesM <- readRDS( file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/GRN/ATACseq_GRNprior/ATACseq_tgenesM.rda")
TF_MeanMed <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFA/TFs_SNSC.podo.MeanMed.rda")

DElist_ATACnets.core <- lapply( 1:4, function(ii){
  datt <- DElist[[ii]]
  TFcore <- TF_MeanMed[ (TF_MeanMed) %in% datt$gName[datt$padj<0.05 & !is.na(datt$padj)]]
  ATACseq_deTF <- ATACseq_tgenesM[ (rownames(ATACseq_tgenesM) %in% (TFcore)) , 
                                     colnames(ATACseq_tgenesM) %in% (TFcore) ]
  # ATACseq_deTF <- ATACseq_deTF[, !colnames(ATACseq_deTF)%in%c("Max", "Nr1d1")]
  # remove TF2 without regulators
  if( sum(colSums(ATACseq_deTF>0)>=1)==1) {
    ATACseq_deTF <- ATACseq_deTF[rowSums(ATACseq_deTF>0)>=1, colSums(ATACseq_deTF>0)>=1 , drop=FALSE]
    
    } else {
      ATACseq_deTF <- ATACseq_deTF[rowSums(ATACseq_deTF>0)>=1, colSums(ATACseq_deTF>0)>=1 ]
      
  }
  # make an edge list from non-square adjacency matrix
  ATACseq_deTF_melt <- reshape2::melt(t(ATACseq_deTF))
  ATACseq_deTF_melt <- ATACseq_deTF_melt[ATACseq_deTF_melt$value>0,]
  
  gg <- graph_from_edgelist( as.matrix(ATACseq_deTF_melt[,-3]) )
  
  V(gg)$size <- sqrt(degree(gg, mode = "out")+1)*5
  return(list(TFcore , gg))
})

### predict TF regulators
{
  XX <- lapply( 1:4 , function(ii){
    datt <- DElist[[ii]]
    DEgenes <- datt$gName[datt$padj<0.05 & !is.na(datt$padj)]
    
    ATACseq_deTF <- ATACseq_tgenesM[ , colnames(ATACseq_tgenesM) %in% (TF_MeanMed) ]
    ATACseq_deTF <- ATACseq_deTF[rowSums(ATACseq_deTF>0)>=1, 
                                 colSums(ATACseq_deTF>0)>=1]
    
    XX <- sapply( seq(ncol(ATACseq_deTF)), function(jj){
      # phyper(Overlap, group2, Total-group2, group1, lower.tail= TRUE)
      group1 <- ATACseq_deTF[,jj]
      group1 <- names(group1)[group1>0]
      Overlap <- length( intersect( DEgenes , group1 ) )
      Total <- nrow(ATACseq_deTF)
      group2 <- length( DEgenes )
      group1 <- length(group1)
      phyper( Overlap, group2 , Total-group2, group1, lower.tail= TRUE)
    })
    names(XX) <- colnames(ATACseq_deTF)
    return(XX)
    })
}

### make one graph for Wt1het.del. and Nphs2 mut.
  {
  TFcore <- Reduce( union, lapply( DElist_ATACnets.core[1:4], `[[` , 1))
  ATACseq_deTF <- ATACseq_tgenesM[ (rownames(ATACseq_tgenesM) %in% (TFcore)) , 
                                   colnames(ATACseq_tgenesM) %in% (TFcore) ]
  # ATACseq_deTF <- ATACseq_deTF[, !colnames(ATACseq_deTF)%in%c("Max", "Nr1d1")]
  # remove TF2 without regulators
  ATACseq_deTF <- ATACseq_deTF[rowSums(ATACseq_deTF>0)>=1, 
                               colSums(ATACseq_deTF>0)>=1]
  
  # make an edge list from non-square adjacency matrix
  ATACseq_deTF_melt <- reshape2::melt(t(ATACseq_deTF))
  ATACseq_deTF_melt <- ATACseq_deTF_melt[ATACseq_deTF_melt$value>0,]
  # make a column which marks edges coming in and out of Wt1
  ATACseq_deTF_melt$value <- ifelse( ATACseq_deTF_melt$Var1=="Wt1" , "out", 
                                     ifelse(ATACseq_deTF_melt$Var2=="Wt1", "in","") )
  ATACseq_deTF_melt$value[ ATACseq_deTF_melt$Var1== "Wt1" & 
                             ATACseq_deTF_melt$Var2=="Wt1" ] <- "in"
  # generate a graph
  gg <- graph_from_data_frame( as.matrix(ATACseq_deTF_melt) )
  # adjust node size
  V(gg)$size <- sqrt(degree(gg)+1)*1.5
  # color nodes, highlight TFs DE in Wt1het.del model and their downstrream TFs.
  XX <- ATACseq_tgenesM[ (rownames(ATACseq_tgenesM) %in% (TF_MeanMed)) , 
                         colnames(ATACseq_tgenesM) %in% union( DElist_ATACnets.core[[2]][[1]],DElist_ATACnets.core[[1]][[1]])  ]
  XX <- XX[rowSums(XX>0)>=1, colSums(XX>0)>=1 ]
  V(gg)$color <- ifelse( V(gg)$name %in% DElist_ATACnets.core[[2]][[1]] , "salmon" ,
                         "palegoldenrod"  )
  # color edges
  E(gg)$color <-  ifelse( E(gg)$value == "out", "#D55E00", 
                          ifelse( E(gg)$value == "in", "#000000","#999999"))
  # plot a graph
  
  plot( gg, edge.arrow.size=.2 ,
        edge.width= 1.2,  
        vertex.label.cex=1.0 ,
        vertex.label.dist=0.1 ,
        # layout = l , 
        # layout=layout_with_mds(gg) , 
        # layout=layout_with_fr(gg, niter=50) ,
        vertex.label = ifelse( degree( gg, mode = "all") > 0, V(gg)$name, NA))
  
  # ggraph::ggraph(DElist_ATACnets.core[[2]][[2]],"stress",bbox = 15) +
  #   geom_edge_link0(edge_colour = "grey66",edge_width = 0.2)+
  #   geom_node_point(aes(fill = exprMat_mod[,i]),shape = 21,size =  sqrt(degree(gg, mode = "all"))*3 )+
  #   geom_node_text(aes(label = ifelse( degree(gg, mode = "all") > 12, V(gg)$name, NA),size = degree(gg)),
  #                  family = "serif",repel = F,
  #                  colour = "darkblue" , size = 8 )+
  #   scale_fill_viridis( limits=c(-2.01,2.01))+
  #   # scale_size(range=c(2,5),guide = FALSE)+
  #   theme_graph()+theme(legend.position = "bottom")
  
}
