#=================================================================#
# ========== FUNCTIONAL ANNOTATION AND PATHWAY ANALYSIS ========= =
#=================================================================#
options(connectionObserver = NULL)
source("/home/tim_nevelsk/PROJECTS/myCode/usefulRfunc.r")
library( biomaRt)
mart_mouse <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "mmusculus_gene_ensembl")
entr2gName <- getBM( attributes=c('external_gene_name', 
                                  "entrezgene_id" ) ,  mart = mart_mouse)
esmbID2gName <- getBM( attributes=c('ensembl_gene_id', 
                                  "external_gene_name" ) ,  mart = mart_mouse)

####======== GO ANALYSIS ============####
### prepare GO annotations
  {
  # #  get all unique genes AND genes of child terms
  # go2entrez <-  as.list(org.Mm.egGO2ALLEGS)
  # go2ensemble <- lapply( go2entrez, function(XX)
  #   unique( entr2gName$ensembl_gene_id[ match( XX , entr2gName$entrezgene_id )] ) )
  # go2gName <- lapply( go2entrez, function(XX)
  #   unique( entr2gName$external_gene_name[ match( XX , entr2gName$entrezgene_id )] ) )
  # # saveRDS(go2ensemble , file="/media/tim_nevelsk/WD_tim/ANNOTATIONS/GO/go2ensemble.rda")
  # saveRDS(go2gName , file="/media/tim_nevelsk/WD_tim/ANNOTATIONS/GO/go2gName.rda")
    go2ensemble <- readRDS( file="/media/tim_nevelsk/WD_tim/ANNOTATIONS/GO/go2ensemble.rda" )
    
    go2gName  <- readRDS( file="/media/tim_nevelsk/WD_tim/ANNOTATIONS/GO/go2gName.rda")
  go2gName <- lapply(go2gName, function(X) X[!is.na(X)])
  go2gName <- go2gName[ sapply(go2gName, length)> 0 ]

  library(GO.db)
  goterms <- Term(GOTERM)
  goont <- Ontology(GOTERM)
  
}

################# GO annotation with Robert function ################ #
  {
  ### prepare gene set by annotating ensembleID with entrezID
    options(connectionObserver = NULL)

  source("/home/tim_nevelsk/PROJECTS/myCode/justGO_ROBERT.R")
  library("org.Mm.eg.db")
    library(GO.db)
    library(DBI)
    
  # # use Robert's function to create sparse binary matrix indicating membership of genes (columns) in GO terms (rows)
  # gomatrix=sf.createGoMatrix()
  # # 31.08.23
  # saveRDS(gomatrix, file="/media/tim_nevelsk/WD_tim/ANNOTATIONS/GO/gomatrix.31.08.23.rda")
  readRDS("/media/tim_nevelsk/WD_tim/ANNOTATIONS/GO/gomatrix.31.08.23.rda")
  # robertGO <- sf.clusterGoByGeneset( gomatrixgomatrix=gomatrix, geneset, universe,
  #                        min.genes=5, cut_max = 10 )

}

####======== GSEA Pathway analysis ####
# prepare Reactome
  {
  # read the file
  reactPath <- fgsea::gmtPathways( "/media/tim_nevelsk/WD_tim/ANNOTATIONS/Pathways/Reactome/ReactomePathways.gmt")
  
  
  humanGen <- unique (Reduce( c , reactPath) )
  mouseGen <- fun_homoTO.FROMmouse(humanGen, TO=T)
  mouseGen$HOMO <- rownames(mouseGen)
  
  reactPath_gName <- lapply( reactPath , function(x){
    X <- as.character( na.omit( unique(mouseGen$MUS[ match( x, mouseGen$HOMO)])))
  } )
  
  # remove untranslated names
  reactPath_gName <- lapply(reactPath_gName, function(X) X <- X[X!="NULL"])
  
  # delete paths with no 3 genes
  reactPath_gName <- reactPath_gName[ unlist( lapply( reactPath_gName, function(x) length(x)>0  ))]
  # remove NAs from gene names 
  reactPath_gName <- lapply(reactPath_gName, function(X) X[!is.na(X)])
  reactPath_gName <- lapply(reactPath_gName, unique)
  
  # create gene set collection
  reactPath_genesets <- lapply( seq(reactPath_gName) , function( jj ) {
    GeneSet( reactPath_gName[[jj]], 
             setName= names(reactPath_gName)[jj] ) } )
  reactPath_genesets <-  GSEABase::GeneSetCollection( reactPath_genesets )
  
  saveRDS( reactPath_genesets , file="reactPath_genesets.13.05.2024.rda")
  
  
}

# prepare KEGG
  {
  ## download directly from KEGG
  library(KEGGgraph)
  kegg_paths <- limma::getGeneKEGGLinks(species.KEGG = "mmu", convert = T)
  X <-  limma::getKEGGPathwayNames(species.KEGG = "mmu", remove.qualifier = T)
  kegg_paths$pNames <- X$Description[ match( sub( "path:","" , kegg_paths$PathwayID) , 
                                             X$PathwayID)]
  kegg_paths$gNames <- entr2gName$external_gene_name[ match(kegg_paths$GeneID , entr2gName$entrezgene_id)]
  # convert long df to a list
  kegg_pathsList <- split(kegg_paths$gNames, kegg_paths$pNames)
  # remove NAs from gene names 
  kegg_pathsList <- lapply(kegg_pathsList, function(X) X[!is.na(X)])
  kegg_pathsList <- lapply(kegg_pathsList, unique)
  
  
  # create gene set collection
  keggPath_genesets <- lapply( seq(kegg_pathsList) , function( jj ) {
    GeneSet( kegg_pathsList[[jj]], 
             setName= names(kegg_pathsList)[jj] ) } )
  keggPath_genesets <-  GSEABase::GeneSetCollection( keggPath_genesets )

  saveRDS( keggPath_genesets , file="keggPath_genesets.13.05.2024.rda")
  
}


####======== Signaling Pathway Impact Analysi (SPIA) ======
#======================================================#
# ############ SPIA via graphite
# ### prepare SPIA
#   {
#   library( "graphite" )
#   library( "SPIA" )
#   library( "org.Mm.eg.db" )
#   options(connectionObserver = NULL)
# 
#   # ## download prepared KEGG and REACTOME pathways to run with graphite
#   path_REACT <-pathways("mmusculus", "reactome")
#   path_REACT <- convertIdentifiers(path_REACT, "ENTREZID")
#   path_KEGG <- pathways("mmusculus", "kegg")
#   # # prepare for SPIA
#   prepareSPIA( path_KEGG, "KEGGAll")
#   prepareSPIA( path_REACT, "reactomeAll")
# 
# 
#   ## databases from 10.2017
#   ## databases from 27-04-2018
#   ## databases from 13-10-2020
#   save( path_KEGG, path_REACT, "/media/tim_nevelsk/WD_tim/ANNOTATIONS/Pathways/SPIA.KeggReactome_13-10-2020.Rdata"  )
#   # load("/media/tim_nevelsk/WD_tim/ANNOTATIONS/Pathways/SPIA.KeggReactome_13-10-2020.Rdata")
#   
#   # ### define SPIA  function
#   # SPIAonList <- function( inputSPIA , ID.key= "SYMBOL" , strict = T , shrunk=T , geneThrsh=0.05 )
#   #   {
#   # 
#   #   ## ID.key  - define what ID type is used:  ENSEMBL or SYMBOL 
#   #   ## inputSPIA - a list of DEseq2 results
#   #   ## strict - a logical argument specifying if all mouse genes (FALSE)
#   #   # or only genes expressed in the experiment (TRUE) are used as the background
#   #   ## shrunk - a logical argument setting if shrunk (TRUE) or non-adjusted (FALSE) lfc valuesare used
#   # 
#   #   require(  "SPIA" )
#   #   library( "graphite" )
#   #   # assign DE results
#   #   x <- inputSPIA
#   #   
#   #   # convert IDs accordingly
#   #   if( ID.key== "ENSEMBL" ){ 
#   #     if(  !exists("geneID2entrez") ){
#   #       geneID2entrez <- biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id"), mart = mart_mouse)
#   #       } 
#   #       # define mapping table
#   #       mappingTab <- geneID2entrez 
#   #       # all genes in experiment
#   #       ALL_genes <- unique( mappingTab$entrezgene[ which(mappingTab$ensembl_gene_id%in%rownames(x) ) ] )
#   #       ALL_genes <- ALL_genes[!is.na(ALL_genes)]
#   #       # modify gene names for compatibility with GRAPHITE: convert to entrez IDs
#   #       x$entrezgene <- mappingTab$entrezgene_id[match( rownames(x), mappingTab$ensembl_gene_id)]
#   #       x <- x[ !is.na(x$entrezgene),]
#   #     } else if( ID.key== "SYMBOL" ) {
#   #       if(  !exists("entr2gName") ){
#   #         entr2gName <- biomaRt::getBM(attributes = c("external_gene_name","entrezgene_id"), mart = mart_mouse)
#   #       } 
#   #       mappingTab <- entr2gName
#   #       # all genes in experiment
#   #       ALL_genes <- unique( mappingTab$entrezgene[ which(mappingTab$external_gene_name%in%rownames(x) ) ] )
#   #       ALL_genes <- ALL_genes[!is.na(ALL_genes)]
#   #       # modify gene names for compatibility with GRAPHITE: convert to entrez IDs
#   #       x$entrezgene <- mappingTab$entrezgene_id[match( rownames(x), mappingTab$external_gene_name)]
#   #       x <- x[ !is.na(x$entrezgene),]
#   #     }
#   #   
#   # 
#   #   if(strict!=T){
#   #     ## define background gene set
#   #     ALL_genes <- union( ALL_genes, unique( mappingTab$entrezgene ) )
#   #     ALL_genes <- paste("ENTREZID", ALL_genes, sep = ":") # modify gene names
#   # 
#   #   } else if (strict==T){
#   #     ## OR define background gene set more strictly as only genes participated in DE analysis
#   #     ALL_genes <- paste("ENTREZID", ALL_genes, sep = ":") # modify gene names
#   #   }
#   # 
#   #   ## select genes DE at geneThrsh significance level
#   #   x <- x[ which( x$padj < geneThrsh), ]
#   #   
#   #   # aggregate LFC for entrez genes with multiple corresponding gene names
#   #   if( isTRUE( shrunk ) ) {
#   #     x <- aggregate( lfc_shrunk~entrezgene,x,sum)
#   #     x <- setNames( x$lfc_shrunk,paste("ENTREZID", x$entrezgene, sep = ":"))
#   #   } else {
#   #     x <- aggregate( log2FoldChange~entrezgene,x,sum)
#   #     x <- setNames( x$log2FoldChange,paste("ENTREZID", x$entrezgene, sep = ":"))
#   #   }
#   #   # # make sure no genes are missing from the background
#   #   # ALL_genes <- union(ALL_genes, names(x) )
#   # 
#   #   ### run SPIA
#   #   SPIA_REACT <- runSPIA(de=x, all=ALL_genes, pathwaySetName="reactomeAll",verbose=F,nB=2000)
#   #   # SPIA_REACT_sig=SPIA_REACT[which(SPIA_REACT$pGFWER<=0.05),]
#   #   SPIA_KEGG <- runSPIA(de=x, all=ALL_genes, pathwaySetName="KEGGAll",verbose=F,nB=2000)
#   #   # SPIA_KEGG_sig=SPIA_KEGG[which(SPIA_KEGG$pGFWER<=0.05),]
#   # 
#   #   return(list(SPIA_KEGG , SPIA_REACT))
#   # }
# 
#   # SPIA_both_NONstrictBkg <- SPIAonList( inputSPIA = , strict = F )
#   # SPIA_both_strictBkg <- SPIAonList(strict = T )
# }

# ### Venn diagram to compare 4 VS 12 weeks pathways
#   { library(VennDiagram)
#   test_compare_gene0.1 <- list(early_FSGS=SPIA_KEGG_w4$Name[which(SPIA_KEGG_w4$pGFWER<0.1)],
#                                late_FSGS=SPIA_KEGG_w12$Name[which(SPIA_KEGG_w12$pGFWER<0.1)])
#   venn.plot <- venn.diagram(test_compare_gene0.1 , NULL, fill=c("darkblue","gray"), alpha=c(0.5,0.5), 
#                             cex = 4, cat.fontface=4, category.names=c("early_FSGS\n N=34", "late_FSGS\n N=26"), main="overlap of KEGG pathways\n identified by SPIA (pGFWER<0.01)")
#   # To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
#   grid.draw(venn.plot)
# }
# 
# ### wrapper function to plot Heatmap results of SPIA 
# path_sel <- read.table(sep = "\t","/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/SPIA/SPIAKEGGREACTOMEOVERLAP_2018.csv")
# path_sel2019 <- path_sel[-c(26,31),]
# # combine REACTOME and KEGG results and add dbID to pathway names
# set4w_both <- rbind(SPIA_REACTlist[[1]],SPIA_KEGGnew[[1]])
# set4w_both$Name <- paste(set4w_both$Name,c(rep("REACTOME", dim(SPIA_REACTlist[[1]])[1]),rep("KEGG",dim(SPIA_KEGGnew[[1]])[1])),sep = "^")
# set12w_both <- rbind(SPIA_REACTlist[[2]],SPIA_KEGGnew[[2]])
# set12w_both$Name <- paste(set12w_both$Name,c(rep("REACTOME", dim(SPIA_REACTlist[[2]])[1]),rep("KEGG",dim(SPIA_KEGGnew[[2]])[1])),sep = "^")
# 
# # load SPIA results 
# ll <- list.files(pattern = "^SPIA_" , "/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/SPIA",full.names = T)
# SPIA_results <- lapply(ll, readRDS)
# SPIA_REACT <-  lapply(SPIA_results, `[[`, 1)
# SPIA_REACT[c(1,2)] <- SPIA_REACT[c(2,1)]
# SPIA_KEGG <-  lapply(SPIA_results, `[[`, 2)
# SPIA_KEGG[c(1,2)] <- SPIA_KEGG[c(2,1)]
# 
# # SPIA_REACT <- list(all_SPIA_REACTallbkg[[1]],all_SPIA_REACTallbkg[[2]],SPIA_REACTnew[[3]],SPIA_REACTnew[[4]])
# # SPIA_KEGG <- list(all_SPIA_KEGGallbkg[[1]],all_SPIA_KEGGallbkg[[2]],SPIA_KEGGnew[[3]],SPIA_KEGGnew[[4]])
# 
# # separate heatmaps 
# # path_sel1 <- path_sel[which(paste(path_sel[,1],path_sel[,2],sep = "^")%in%union_tab$union[which(union_tab[,2]<=0.1 & union_tab[,4]<=0.1)]),]
# # path_sel2 <- path_sel[which(paste(path_sel[,1],path_sel[,2],sep = "^")%in%union_tab$union[which(union_tab[,2]<0.1 & union_tab[,4]>0.1)]),]
# # path_sel3 <- path_sel[which(paste(path_sel[,1],path_sel[,2],sep = "^")%in%union_tab$union[which(union_tab[,2]>0.1 & union_tab[,4]<0.1)]),]
# 
# nameCompare <-  c("FSGS_4w","FSGS_12w","aging","dis.prog.","dis.prog.(-aging)")
# 
# plotSPIA <- function( set1, set2, SPIA_KEGGlist=NULL, SPIA_REACTlist=NULL, nameL = 50 , DendroBar=F, ringPlot=F,  scaleR=c(-10,10),
#                       alpha=0.1, path_type="path_type", selected=NULL, C_height = 5 , orderLFC=F , colab=nameCompare)
#  { 
#   require(ggplot2)
#   require(reshape)
#   require(plyr)
#   
#   if (is.null(SPIA_KEGGlist)){
# 
#    if (is.null(selected)){
#     set1_sig <- set1[which(set1$pGFWER<alpha),]
#     set2_sig <- set2[which(set2$pGFWER<alpha),]
#   } else {
#     selected <- paste(selected[,1],selected[,2],sep = "^")
#     names(set1) <- names(set2) 
#     set1_sig <- set1[which(set1[,1]%in%selected),] 
#     set2_sig <- set2[which(set2[,1]%in%selected),] 
# 
#   }
#   
#   union <- union(set1_sig$Name, set2_sig$Name)
#   
#   union_tab <- cbind(union,set1[(match(union,set1$Name)),c(9,10)],set2[(match(union,set2$Name)),c(9,10)])
#   union_tab[is.na(union_tab[,2]),2] <- 1  # eleminate possible NAs
#   union_tab[is.na(union_tab[,4]),4] <- 1
#   union_tab[is.na(union_tab[,3]),3] <- "Inhibited"
#   union_tab[is.na(union_tab[,5]),5] <- "Inhibited"
#   names(union_tab)=c(path_type,"early_FSGS_pGFWER","set1_status","late_FSGS_pGFWER","set2_status")
#   
#   X <- union_tab
#   names <- gsub("_pGFWER","",names(X)[c(1,2,4)])
#   
#   x <- revalue(as.factor(X[,3]), c("Activated"=-1, "Inhibited"=1))
#   x =as.numeric(levels(x))[x]
#   y <- revalue(as.factor(X[,5]), c("Activated"=-1, "Inhibited"=1))
#   y= as.numeric(levels(y))[y]
#   
#   inter_toPlot=data.frame(both=X[,1],FSGS_4w=log10(X[,2])*x, FSGS_12w=log10(X[,4])*y)
#   } else { 
#     if(is.null(selected)) {stop("provide a list of selected patways to plot SPIA results for all comparisons")
#       } else{
#       # ll <- length(SPIA_KEGGlist[[1]][match( selected[,1] , SPIA_KEGGlist[[1]]$Name),1] )
#       path_union <-   paste(selected[,1],selected[,2],sep = "^")
#       
#       # check if length of Reactome and SPIA are equal if no throw an error
#       if(length(SPIA_KEGGlist)==length(SPIA_REACTlist)) {
#         
#       setAll_sig <- lapply(seq(SPIA_KEGGlist) , function(x){
#           set_sig <- rbind( SPIA_KEGGlist[[x]] , SPIA_REACTlist[[x]] )
#           set_sig$Name <- c( paste( set_sig$Name[1:dim(SPIA_KEGGlist[[x]])[1]], '^KEGG' , sep = "") , paste(set_sig$Name[(dim(SPIA_KEGGlist[[x]])[1]+1):dim(set_sig)[1]], '^REACTOME' , sep = "") )
#           set_sig <- set_sig[ match( path_union , set_sig$Name),] 
#           return(set_sig)
#         })
#         
# 
#         ### union of all comparisons
#         union_tab <- lapply(setAll_sig , function(x) x[(match(path_union,x[,1])),c(9,10)])
#         union_tab <- Reduce( cbind , union_tab) 
#         union_tab <- cbind.data.frame( path_union, union_tab )
#         
#         union_tab[ , seq( 2 ,(dim(union_tab)[2]-1), by=2)][is.na(union_tab[ , seq( 2 ,(dim(union_tab)[2]-1), by=2)] )] <- 1  # eleminate possible NAs
#         union_tab[ , seq( 3 ,(dim(union_tab)[2]), by=2)][is.na(union_tab[ , seq( 3 ,(dim(union_tab)[2]), by=2) ] )] <- "Inhibited"
#         
#         X <- union_tab
#         # revalue in a list
#         X[ , seq( 3 , (dim(X)[2]) , by=2) ] <- lapply( X[, seq( 3 , (dim(X)[2]) , by=2)] , function(y) {
#           x <- revalue(as.factor(y), c("Activated"=-1, "Inhibited"=1)) 
#           x <- as.numeric(levels(x))[x]
#           return(x) })
#         
#         inter_toPlot <- cbind.data.frame( X[,1] , log10(X[ , seq( 2 ,(dim(union_tab)[2]-1), by=2)] ) * X[, seq( 3 ,(dim(union_tab)[2]), by=2)] )
#         names(inter_toPlot) <- c( "both" , colab )
#       } else stop("lists of SPIA results for KEGG and REACTOME should be of the same length")
#      
#     }
#     
#   }
# 
#     # round significance values (in case we want to display them)
#     inter_toPlot[,-1] <- round((inter_toPlot[,-1]),2)
#     inter_toPlot[,1] <- sapply(inter_toPlot[,1], function(x) wrap_text(as.character(x),n=nameL) )
#     
#     if (!is.null(scaleR)) {
#       if (length(scaleR)==2 & class(scaleR)=="numeric"){
#         #inter_toPlot <- cbind.data.frame(both=inter_toPlot[,1],t(scale(t(inter_toPlot[,-1])) ))
#         inter_toPlot[,-1][inter_toPlot[,-1] < (scaleR[1])] <- scaleR[1]
#         inter_toPlot[,-1][inter_toPlot[,-1] > (scaleR[2])] <- scaleR[2]
#       } else stop("scaleR argument is incorrect: set to NULL or supply a numeric vector of length 2 with lower and upper limits for lfc scale")
#         
#   }
#    
#     if(orderLFC==F){
#       ordR <- hclust( dist(inter_toPlot[,-1] , method = "euclidean"), method = "ward.D" )$order
#     } else  ordR <- order(rowMeans(inter_toPlot[,-1]))
#     
# 
#     # cluster pathways
#     # melt for use with ggplots
#     melted_cormat <- melt(inter_toPlot, na.rm = TRUE)
#     melted_cormat[,1] <- factor(melted_cormat[,1], levels = unique(melted_cormat[,1]))
#     melted_cormat[,4] <- abs( melted_cormat[,3])
# 
#     # order rows
#     melted_cormat[,1] <- factor( melted_cormat[,1] , levels = inter_toPlot[,1][ordR], labels =   inter_toPlot[,1][ordR] )
# 
#     if(ringPlot==F) {
#       
#       # if(value==T) # use of values in tiles is depricated for now due to no utility
#       #   p1 <- ggplot(data = melted_cormat, aes(variable, melted_cormat[,1], fill = value)) + geom_tile(color = "white",stat = "identity",size=0)+
#       #   geom_text(aes(label=V4),cex = 2.5) +  scale_y_discrete(expand=c(0,0) , position = "right") + # to show values in tiles
#       #   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#       #                        midpoint = 0, limit = c(min(melted_cormat$value),max(melted_cormat$value),labels=NULL), space = "Lab", 
#       #                        name="Direction of pathway change\n pseudo-colours") +
#       #   theme_minimal()+  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle = 45,hjust=1, size=14),
#       #                           axis.title.y=element_blank())  + coord_fixed(ratio=0.2)
#       # else
#       
#       p1 <- ggplot(data = melted_cormat, aes(variable, melted_cormat[,1], fill = value)) + geom_tile(color = "white",stat = "identity",size=0)+
#         scale_y_discrete(expand=c(0,0) , position = "right") + scale_x_discrete(expand=c(0,0) ) +
#         scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                              midpoint = 0, labels=abs , space = "Lab", 
#                              limit = c(min(melted_cormat$value),max(melted_cormat$value)),
#                              name="Pathway change\nlog10(pGFWER)\nin pseudo-colours") +
#         theme_minimal()+  theme(axis.title=element_blank(), 
#                                 element_line(colour=NULL) ,
#                                 axis.ticks= element_blank(),
#                                 axis.text.x=element_text(angle = 90, hjust=1 , vjust=0.5 , size=14 ),
#                                 axis.text.y = element_text(size=10 ), axis.title.y=element_blank())  + coord_fixed(ratio=0.2)
#       
#       if(DendroBar==T) {
#         
#         # make cluster bar in addition to a heatmap
#         require(ggdendro)
#         require(gridExtra)
#         require(ggpubr)
#         
#       hc <- hclust( dist(as.data.frame( inter_toPlot[,-1] , row.names = as.character( inter_toPlot[,1]) ) , method = "euclidean"), method = "ward.D" )
#       df2 <- data.frame( cluster=cutree(hc,C_height) , states=factor(hc$labels,levels=hc$labels[hc$order]))
#       # plot
#       p2 <- ggplot( df2, aes(x=1,y=states,fill=factor(cluster))) + geom_tile()+
#         scale_y_discrete(expand=c(0,0))+
#         scale_x_discrete(expand=c(0,0),name = NULL ) +
#         theme( axis.title.y = element_blank() , element_line(colour="white") ,
#                axis.ticks= element_blank(),
#                axis.text= element_blank(),
#                axis.title.x = element_text(face="bold" , angle = 90 , size=14 ,  vjust = 0.5),
#                legend.position="none" )
#       ## arrange grids
#       ggpubr::ggarrange( p2, p1, ncol=2 , widths=c(1,20) , align = "h")
#       
#      } else(p1)
#     } else {
#       require(forcats)
#       sequence_length = length(unique(melted_cormat[,1]))
#       first_sequence = c(1:(sequence_length%/%2)) 
#       second_sequence = c((sequence_length%/%2+1):sequence_length) 
#       first_angles = c(90 - 180/length(first_sequence) * first_sequence)
#       second_angles = c(-90 - 180/length(second_sequence) * second_sequence)
#       
#       melted_cormat$variable <- fct_rev(melted_cormat$variable)
#         
#       ggplot2::ggplot(data = melted_cormat, aes(x=variable, y=melted_cormat[,1], fill = value )) + geom_tile(color = "white",stat = "identity",size=0)+
#         scale_x_discrete(expand=c(1,1)) +  
#         scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                              midpoint = 0, labels=abs , space = "Lab", 
#                              #limit = c(min(melted_cormat$value),max(melted_cormat$value)),
#                              name="Pathway change\nlog10(pGFWER)\nin pseudo-colours") + 
#          theme_minimal()+  
#           theme(axis.title=element_blank(), plot.margin=unit(c(5,5,5,5),"cm") ,
#                 panel.border = element_blank(),
#                 #panel.grid.major = element_blank(),
#                                 element_line(colour=NULL) ,
#                                 axis.ticks= element_blank(),
#                                 axis.text.x = element_text(angle = c(first_angles,second_angles), hjust=rep(1,length(melted_cormat[,1])) , vjust=-1 , size=10 ),
#                                 axis.text.y = element_blank(), axis.title.y=element_blank()) +  
#         coord_polar(theta="y", clip="off") # circularize plot, dont clip pathway names: clip="off"
#         
#         }
# }
# 
# 
# # makle a plot (circular)
# plotSPIA(   set1= set4w_both, set2=set12w_both,  DendroBar=F, ringPlot=F, scaleR=NULL, path_type="path_type", 
#             alpha=0.1 , C_height = 5, nameL = 200 ,  orderLFC=F )
# # non circular
# plotSPIA(   SPIA_KEGGlist=SPIA_KEGG , SPIA_REACTlist=SPIA_REACT ,  DendroBar=F, ringPlot=F, scaleR=c(-15,15) , path_type="path_type", 
#             selected = path_sel2019, C_height = 5, nameL = 200 ,  orderLFC=F )
# 
# ### heatmap for selected paths
#  {    
# # SPIA_4w <- rbind(SPIA_KEGG_w4,SPIA_REACT_w4)
# #   SPIA_4w[,11] <- c(rep("KEGG",dim(SPIA_KEGG_w4)[1]),rep("REACT",dim(SPIA_REACT_w4)[1]))
# #   SPIA_4w <- subset(SPIA_4w, (SPIA_4w$Name%in%path_sel$V1[which(path_sel$V2=="KEGG")] & SPIA_4w$V11=="KEGG") | (SPIA_4w$Name%in%path_sel$V1[which(path_sel$V2=="Reactome")] & SPIA_4w$V11=="REACT"))
# #   SPIA_4w$Name[SPIA_4w$V11=="REACT"] <- paste(SPIA_4w$Name[SPIA_4w$V11=="REACT"],"*",sep = "")
# #   SPIA_12w <- rbind(SPIA_KEGG_w12,SPIA_REACT_w12)
# #   SPIA_12w[,11] <- c(rep("KEGG",dim(SPIA_KEGG_w12)[1]),rep("REACT",dim(SPIA_REACT_w12)[1]))
# #   SPIA_12w <- subset(SPIA_12w, (SPIA_12w$Name%in%path_sel$V1[which(path_sel$V2=="KEGG")] & SPIA_12w$V11=="KEGG") | (SPIA_12w$Name%in%path_sel$V1[which(path_sel$V2=="Reactome")] & SPIA_12w$V11=="REACT"))
# #   SPIA_12w$Name[SPIA_12w$V11=="REACT"] <- paste(SPIA_12w$Name[SPIA_12w$V11=="REACT"],"*",sep = "")
# # 
# #   # prepare data for plotting
# #   library(plyr)
# #    union <- union(SPIA_4w$Name,SPIA_12w$Name)
# #   union_tab <- cbind(union,SPIA_4w[(match(union,SPIA_4w$Name)),c(9,10)],SPIA_12w[(match(union,SPIA_12w$Name)),c(9,10)])
# #   union_tab[is.na(union_tab[,2]),2] <- 1  # eleminate possible NAs
# #   union_tab[is.na(union_tab[,4]),4] <- 1
# #   union_tab[is.na(union_tab[,c(3)]),c(3)] <- "Inhibited"
# #   union_tab[is.na(union_tab[,c(5)]),c(5)] <- "Inhibited"
# #   names(union_tab)=c("pathway","early_FSGS_pGFWER","set1_status","late_FSGS_pGFWER","set2_status")
# #   union_tab[union_tab$early_FSGS_pGFWER<0.1 & union_tab$late_FSGS_pGFWER<0.1,6] <- "both"
# #   union_tab[union_tab$early_FSGS_pGFWER<0.1 & union_tab$late_FSGS_pGFWER>=0.1,6] <- "early_only"
# #   union_tab[union_tab$early_FSGS_pGFWER>=0.1 & union_tab$late_FSGS_pGFWER<0.1,6] <- "late_only"
# #   
# #   #  x <- union_tab[union_tab$V6=="both",]
# #   # names <- gsub("_pGFWER","",names(x)[c(1,2,4)])
# #   # inter_toPlot=x[order(x[,3],x[,2]),]     # order according to status and p-value of the first condition
# #   # x=revalue(as.factor(inter_toPlot[,3]), c("Activated"=-1, "Inhibited"=1))
# #   # x=as.numeric(levels(x))[x]
# #   # y=revalue(as.factor(inter_toPlot[,5]), c("Activated"=-1, "Inhibited"=1))
# #   # y=as.numeric(levels(y))[y]
# #   # inter_toPlot=data.frame(inter_toPlot[,1],log10(inter_toPlot[,2])*x,log10(inter_toPlot[,4])*y)
# #   # names(inter_toPlot) <- names
# #   # inter_toPlot=inter_toPlot[order(inter_toPlot[,2],decreasing = T),]
# #   # inter_toPlot[,2:3] <- round((inter_toPlot[,2:3]),2)
# #   
# #   union_tab_trans <- 
# #   
# #   melted_cormat <- melt( union_tab[,c(1,2,4)] , na.rm = TRUE)
# #   # melted_cormat[,1] <- factor(melted_cormat[,1], levels = unique(melted_cormat[,1]))
# #   # melted_cormat[,4] <- abs(melted_cormat[,3])
# #   melted_cormat$variable <- factor(melted_cormat$variable, labels = c("early FSGS","late FSGS"))
# #   
# #   
# #   labl <- c(20,"10 down-regulated",5,0,5,"10 up-regulated")
# #   
# #   ordR <- hclust( dist(XX , method = "euclidean"), method = "ward.D" )$order
# #   ordC <- hclust( dist( t(XX) , method = "euclidean"), method = "ward.D" )$order
# #   
# #   XX.m <- melt( as.matrix( XX ) )
# #   # # add significance marks to tiles
# #   XX.m$padj <- melt(HomoMouse_qval)$value
# #   XX.m$padj <-  ifelse( is.na(XX.m$padj) , "", ifelse( XX.m$padj > 0.1 , "" , ifelse( XX.m$padj > 0.01 ,  "*" , "**") ) )
# #   
# #   # order rows
# #   XX.m$X1 <- factor( XX.m$X1, levels = rownames(XX)[ordR], labels =  rownames(XX)[ordR] )
# #   # order columns
# #   XX.m$X2 <- factor( XX.m$X2, levels = colnames(XX)[ordC], labels =  colnames(XX)[ordC] )
# #   
# #   g1<- ggplot(data = melted_cormat, aes(variable, melted_cormat[,1], fill = value)) + geom_tile(color = "white",stat = "identity",size=0)+
# #     scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
# #                          midpoint = 0, limit = c(min(melted_cormat$value),max(melted_cormat$value)),labels=labl, space = "Lab", 
# #                          name="-log10*q-value\n of pathway change\n in pseudo-colours") +
# #     theme_minimal()+  theme(text = element_text(size=14),axis.title.x=element_blank(),axis.text.x=element_text(angle = 45,hjust=1),
# #                             axis.title.y=element_blank())  + coord_fixed(ratio=0.2)
# #   # p1 <- as.grob(g1)
# #   # p2 <-  as.grob(g2)
# #   # p3 <-  as.grob(g3)
# #   # grid.arrange(grobs=list(p1, p2, p3),ncol=1,heights = c(2, 1, 1))
#  } 
#    
# ### venn of pathway overlap
#   library(VennDiagram)
#     test_compare_gene0.1 <- list(early_FSGS=SPIA_4w$Name[which(SPIA_4w$pGFWER<0.1)],
#                                  late_FSGS=SPIA_12w$Name[which(SPIA_12w$pGFWER<0.1)])
#     venn.plot <- venn.diagram(test_compare_gene0.1 , NULL, fill=c("green","salmon"), alpha=c(0.5,0.5),
#                               cex =4,cat.fontface=4, cat.cex=2,cat.just=list(c(0.2,-0.35),c(0.8,-0.35)),category.names=c("early FSGS\n N=26", "late FSGS\n N=25"))
#     # To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
#     grid.draw(venn.plot)
# 
# 
# ### pathview visualization of KEGG pathways highlited by SPIA ###
#   { 
#     
#     library(plyr)
#   # make an LFC table with ENTREZ ids
#   LFC_table <- data.frame(WT1_4w_LFC=WT1_4w$log2FoldChange,WT1_12w_LFC=WT1_12w$log2FoldChange,row.names = rownames(WT1_4w))
#   x<- merge(LFC_table,geneID2entrez,by.x = 0, by.y = "ensembl_gene_id",all = F)
#   x <- x[!is.na(x$entrezgene),]
#   x <- ddply(x[,-1],"entrezgene",numcolwise(sum))
#   LFC_tableKEGG <- data.frame(WT1_4w_LFC=x$WT1_4w_LFC,WT1_12w_LFC=x$WT1_12w_LFC,row.names = x$entrezgene)
#   # select KEGG pathways to vizualize
#   set1_sig <- set1[which(set1$pGFWER<0.1),]
#   set2_sig <- set2[which(set2$pGFWER<0.1),]
#   union <- union(set1_sig$Name,set2_sig$Name)
#   
#   setwd("/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/SPIA/KEGG_PATHVIEW/")  
#   for(i in 1:length(union)){
#     kegg_path <- path_KEGG[[which(names(path_KEGG)==union[i])]]
#     
#     # pathview(gene.data=LFC_tableKEGG, pathway.id=sub("mmu:","",kegg_path@id), species="mmu",kegg.native=T, out.suffix = paste(kegg_path@title,"FSGS",sep="."),
#     #          kegg.dir = "/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/RNAseq/SPIA/PATHVIEW/KEGG_path",bins=list(gene=20, cpd=10))
#     pathview(gene.data=LFC_tableKEGG, pathway.id=sub("mmu:","",kegg_path@id), species="mmu",kegg.native=F, out.suffix = paste(kegg_path@title,"FSGS",sep="."),
#              kegg.dir = "/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/RNAseq/SPIA/KEGG_PATHVIEW/KEGG_path",bins=list(gene=20, cpd=10))
#   }
# }  
# 
### 2-way evidence SPIA plot
      spia2way_ggplot <- function( pathSPIA , namesSPIA , pathtype="REACTOME",
                                   SPIA_list = SPIA_REACT,
                                   combinemethod = "fisher", threshold = 0.1)
        {
        library(SPIA)
        library(ggplot2)
        library(ggrepel)
        # helper function from SPIA
        getP2<-function( pG, combine=combinemethod )
        {
          #given a pG returns two equal p-values such as   combfunc(p1,p2)=pG
          if(combine=="fisher"){
            ch=qchisq(pG,4,lower.tail = FALSE)
            return(sqrt(exp(-ch/2)))
          }

          if(combine=="norminv"){
            return(pnorm(qnorm(pG)*sqrt(2)/2))

          }
        }

        PPlist <- list()

        PPlist <- lapply( seq(SPIA_list), function(i)
        {

          x <- SPIA_list[[i]]
          x$pPERT[is.na(x$pPERT)] <- 1
          x$pNDE[is.na(x$pNDE)] <- 1

          # check if there are any terms with padj < threshold
          # if none - return NA instead of a plot
          if( nrow(x[x$pGFdr<=threshold,])==0 ) return(NA) else{
            #determine what combine method was used to convert ph and pb into pG
            # combinemethod=ifelse(sum(combfunc(x$pPERT,x$pNDE,"fisher")==x$pG)>sum(combfunc(x$pPERT,x$pNDE,"norminv")==x$pG),"fisher","norminv")

            tr1<-threshold/dim(na.omit(x))[1]
            trold=tr1
            tr2<-max(x[,"pG"][x[,"pGFdr"]<=threshold])
            if(tr2<=trold) tr2=trold*1.03
            # color points based on significance
            # oks<-x[,"pGFWER"]<=threshold
            # oks2<-x[,"pGFdr"]<=threshold
            x$colourSig <- ifelse( x$pG < tr1, "red",  ifelse( x$pG < tr2, "blue","dark grey") )
            x$colourDir <- ifelse( (x$Status =="Activated" & x$pG < tr2 ), "red", ifelse( x$pG < tr2, "blue","dark grey") )

            x$NameSig <-  ifelse( x$pG < tr1, x$Name, "" )

            # x$pSizeSQRT <- sqrt(x$pSize)

            # set hard limits on X and Y coordinates
            # x$pPERT <- ifelse( -log(x$pPERT) > 14.9 , 0.00000001 ,  x$pPERT )
            # x$pNDE <- ifelse( -log(x$pNDE) > 14.9 , 0.00000035 ,  x$pNDE )
            PP <-  ggplot(x, aes(x=-log(pNDE),y=-log(pPERT)) )+ geom_point( aes(x =-log(pNDE), y=-log(pPERT),colour=colourDir , size = pSize ) , alpha=0.8 ) +
              scale_size( range = c(1,20)) +
              scale_colour_manual(name = 'Pathway status', values = setNames(c("dark grey","blue","red") ,c("dark grey","blue","red")), labels=c('non-significant','down-regulated','up-regulated'))+
              geom_abline(aes(intercept = -log(getP2(tr1,combinemethod)^2), slope = -1), color = "green" ,lwd=1)+
              geom_abline(aes(intercept = -log(getP2(tr2,combinemethod)^2), slope = -1), color = "grey50" ,lwd=1) +
              # ylim(0, 15) + xlim(0, 15) +
              ggtitle(paste(names(SPIA_list)[i], pathtype, sep = " "))+

              theme_minimal() + theme(plot.title = element_text(face="bold",size=20),
                                      axis.text=element_text(size=20),
                                      axis.title=element_text(size=20,face="bold") ,
                                      legend.text = element_text(size = 20),
                                      legend.title = element_text(size=20),
                                      legend.key.size = unit(2, "cm"),
                                      strip.text.y = element_text(size = 20))+
              guides(colour = guide_legend(override.aes = list(size=20))) + # control size of of points in legend
              theme(legend.position = "none") +
              geom_text_repel( size = 6 , point.padding =1 ,
                               aes(label = x$NameSig),box.padding   = 1,
                               segment.color = 'grey50', max.overlaps =100 )  #ensure that labels are not jammed together
            # print(PP)
            # pdf(file = paste(pathSPIA,namesSPIA[i],"_SPIAevid_",pathtype,".pdf",sep = ""), width = 15, height = 12)
            # print (PP)
            # dev.off()
            #
            # png(file = paste(pathSPIA,namesSPIA[i],"_SPIAevid_",pathtype,".png",sep = ""), width = 1000, height = 700)
            # print (PP)
            # dev.off()
            return(PP)
          }



        })
        return(PPlist)
      }


# ### compare 2 conditions using SPIA
#     {
#       
#       SPIA_Nphs2 <- readRDS("/data/user/tpadvits/PROJECTS/MARTIN_PJ/RNAseq_analysis/NPHS2/SPIA/Nphs2_SPIA.nonstrict.rda")
#       # path_sel <-  read.table(sep = "\t","SPIAKEGGREACTOMEOVERLAP_2018.csv")
#       # path_sel <- paste(path_sel[,1],path_sel[,2],sep = "^")
#       
#     # define string wrapper function
#       wrap_text <- function(string, n=40) 
#       { require(stringr)
#         if(nchar(string) > n){
#           spaces <- str_locate_all(string, " ")[[1]][,1]
#           chars  <- nchar(string)
#           for(i in 1:floor(chars/n)) {
#             s <- spaces[which.min(abs(spaces - n*i))]
#             substring(string, s, s) <- "\n "
#           }
#           return(string)
#         } else string <- string
#         
#       }
#       
#       plotSPIA_2D <- function( SPIA_results=SPIA_Nphs2 , nameL = 50 ,
#                             alpha=0.05 , selected=NULL , Nsiglabel=20 , plotTitle="")
#       {    
#         #path_sel <- path_sel[-c(26,31),]
#         set4w_both <- rbind( SPIA_results[[2]][[1]] , SPIA_results[[1]][[1]] ) 
#         set4w_both$Name <- paste( set4w_both$Name, c(rep("REACTOME", dim( SPIA_results[[2]][[1]])[1] ),rep("KEGG",dim(SPIA_results[[1]][[1]])[1] )),sep = "^") 
#         set12w_both <- rbind( SPIA_results[[2]][[2]] , SPIA_results[[1]][[2]] ) 
#         set12w_both$Name <- paste(set12w_both$Name,c(rep("REACTOME", dim( SPIA_results[[2]][[2]])[1] ),rep("KEGG",dim(SPIA_results[[1]][[2]])[1] )),sep = "^")
#         
#       set1_sig <- set4w_both[which(set4w_both$pGFWER<alpha),]
#       set2_sig <- set12w_both[which(set12w_both$pGFWER<alpha),]
#       
#       union <- union(set1_sig$Name, set2_sig$Name)
#       
#       union_tab <- cbind( union, set4w_both[(match(union,set4w_both$Name)),c(9,10)],set12w_both[(match(union,set12w_both$Name)),c(9,10)])
#       union_tab[is.na(union_tab[,2]),2] <- 1  # eleminate possible NAs
#       union_tab[is.na(union_tab[,4]),4] <- 1
#       union_tab[is.na(union_tab[,3]),3] <- "Inhibited"
#       union_tab[is.na(union_tab[,5]),5] <- "Inhibited"
#       
#       X <- union_tab
#       x <- revalue(as.factor(X[,3]), c("Activated"=-1, "Inhibited"=1))
#       x <- as.numeric(levels(x))[x]
#       y <- revalue(as.factor(X[,5]), c("Activated"=-1, "Inhibited"=1))
#       y <- as.numeric(levels(y))[y]
#       
#       inter_toPlot=data.frame( both=X[,1],FSGS_4w=log10(X[,2])*x, FSGS_12w=log10(X[,4])*y)
#       inter_toPlot$shapeD <- as.factor( ifelse( abs(inter_toPlot$FSGS_4w) > 1 & abs(inter_toPlot$FSGS_12w) > 1, 'change in both' ,ifelse(  abs(inter_toPlot$FSGS_4w) < 1 & abs(inter_toPlot$FSGS_12w) > 1 , '12w only' , '4w only' ) ) )
#       
# 
#       inter_toPlot$colorD <- factor(  inter_toPlot$shapeD , levels = c('change in both','4w only', '12w only') , labels = c('green','magenta','orange'))
#    
#       
#       inter_toPlot$both <- ifelse( inter_toPlot$both%in% union( set1_sig$Name[order(set1_sig$pGFWER)][ seq(Nsiglabel) ] , set2_sig$Name[order(set2_sig$pGFWER)][ seq(Nsiglabel) ]) , as.character(inter_toPlot$both) , "")
#       
#       # wrap long pathway names
#       inter_toPlot[,1] <- sapply(inter_toPlot[,1], function(x) wrap_text(as.character(x),n=nameL) )
#       
#       
#       ggplot(data = inter_toPlot, aes(x=FSGS_4w, y=FSGS_12w, shape=shapeD)) + geom_hline( yintercept=0 ) + geom_vline( xintercept=0 )  +
#         geom_point(size=5,  fill=inter_toPlot$colorD  ) + 
#         scale_shape_manual("", values = c( 'change in both'=21,'4w only'=22, '12w only'=24)) + ggtitle(plotTitle) +
#         geom_text_repel(data=inter_toPlot,aes(label=both), box.padding   = 0.35, point.padding = 0, segment.color = 'grey50', size=4)  + #ensure that labels are not jammed together
#         theme_minimal() +  theme(axis.text=element_text(size=20), plot.title = element_text(size = 24, face = "bold", hjust = 0.5) ,
#                                  axis.title=element_text( size=24,face="bold"), legend.text = element_text(size = 24),legend.title = element_text(size=24) )+
#         scale_size(range = c(5, 6)) + theme(legend.key.size = unit(2.5, "cm")) + labs(colour = "condition") # +
#       }
#       
#       
#     }
#     
####======== ENRICHr - web-based tool ===================
## http://amp.pharm.mssm.edu/Enrichr/



####======== barplots for annotation results ==========
 # barplot for robertGO results     
      GOrobert_barplot <-  function(iid=iids, 
                                    GOtoPlot=topN_GOrb,
                                    datGO = GOrobert_res ,
                                    datLFC= GO_meanLFC,
                                    datOrder,
                                    labelTRUNK=T,
                                    typeSig=T)
      {
        ## * iids - a numeric vector, defining which elements of datGO list to plot as columns of the barplot
        ## * GOtoPlot - a charchter vector of GO IDs for which datLFC will be plotted 
        ## * datGO - a list, where each element contains results of Robert's sf.createGoMatrix() funcion
        ## * datLFC - a matrix of mean LFCs (or other measures of GO changes) 
        ## that will be used as a color of bars, if not provided log2FC of enrichment will be used
        ## * oorder - a charachter vector, defining order in which datasets (in columns) appear on the plot
        ## * labelTRUNK - logical value, if labels should be shorten
        ## * typeSig - logical value, whether to split rows by Ontology type of universality of terms
        require(colorspace)

        datTOplot <- Reduce( rbind , lapply( iid , function(ii){
          print(ii)
          datt <- datGO[[ii]]
          # select results for a given GO term, even if it's not a primary term in a given comparison
          datt <- datt$results[ match( GOtoPlot , datt$results$GO.ID )  , ]
          datt$dataset <- names(datGO)[ii]
          
          print( dim(datt) )
          
          ### add LFC or use log2Enrichment of enrichment
          if( !exists("datLFC") ) {
            datt$LFC <- datt[GOtoPlot , 'log2Enrichment']
            } else datt$LFC <- datLFC[GOtoPlot, ii]
          

          # make sore no terms are missing 
          datt$Ontology <- Ontology(GOTERM)[GOtoPlot]
          datt$Term <- Term(GOTERM)[GOtoPlot]
          datt$GO.ID <- GOtoPlot
          datt$Fisher[ is.na(datt$Fisher)] <-  1
          if( isTRUE(typeSig) ) datt$Ontology <- names(GOtoPlot)
          return(datt)
        }))
        
        
        # reorder levels
        datTOplot$dataset <- factor( datTOplot$dataset,
                                     levels = datOrder)
        # truncate labels
        if( isTRUE(labelTRUNK)) datTOplot$Term <- sapply( datTOplot$Term , str_trunc, width = 40)
        
        gg<- ggplot( datTOplot, aes(x=-log10(Fisher), 
                                    y=reorder( Term, -log10(Fisher) ,FUN = mean, na.rm=T ),
                                    fill=LFC ) ) + 
          scale_fill_continuous_divergingx(palette = "Spectral",
                                           na.value="grey",
                                           rev =T, mid = 0 ) +  
           scale_size(range = c(0, 15))+ ylab("GO terms") +
          # geom_vline(xintercept=0, linetype=2,lwd=2)+ # add vertical line for lfc=0
          geom_bar(stat = "identity") + 
          facet_grid(cols = vars(dataset), rows = vars(Ontology), 
                     scales = "free_y", space = "free")+
          theme_bw() + theme( text = element_text(size = 20 ),
                              axis.text.x = element_text(size = 16))
        
        return(gg)
      }
    
    # clustered barplot for robertGO results     
      GOrobert_barplot_clustered <-  function(iid=iids, 
                                    GOtoPlot=topN_GOrb,
                                    datGO = GOrobert_res ,
                                    datLFC= GO_meanLFC,
                                    datOrder,
                                    labelTRUNK=T,
                                    typeSig=T)
        {
        ## * iids - a numeric vector, defining which elements of datGO list to plot as columns of the barplot
        ## * GOtoPlot - a charchter vector of GO IDs for which datLFC will be plotted 
        ## * datGO - a list, where each element contains results of Robert's sf.createGoMatrix() funcion
        ## * datLFC - a matrix of mean LFCs (or other measures of GO changes) 
        ## that will be used as a color of bars, if not provided log2FC of enrichment will be used
        ## * oorder - a charachter vector, defining order in which datasets (in columns) appear on the plot
        ## * labelTRUNK - logical value, if labels should be shorten
        ## * typeSig - logical value, whether to split rows by Ontology type of universality of terms
        require(colorspace)
        library(GeneOverlap)
        
        ### cluster for similarity of paths
        {
          # table of similarity
          GOlist <-go2gName[names(go2gName) %in% GOtoPlot]
          MM<- newGOM( GOlist , GOlist ,genome.size=10000)
          MMpval <- getMatrix(MM, name="pval")
          MMpval[MMpval==0] <- NA
          diag(MMpval) <- NA
          MMpval <- -log10(MMpval)
          
          ## plot
          toPlot <- MMpval  

          # pdf( width =  12 , height = 9, file="gseaKEGG.REACT.top3_pvalHeat.clust.pdf")
          # a<-dev.cur()
          # png( width =  1200 , height = 800 , file="gseaKEGG.REACT.top3_pvalHeat.clust.png")
          # dev.control("enable")
          gsea.order <- gplots::heatmap.2( toPlot,
                                           col = RColorBrewer::brewer.pal( name = "Greens", n=9), 
                                           trace = "none", symkey=T,
                                           na.color = "#00441B",
                                           dendrogram = "row",margins = c(2,20) ,
                                           labCol = FALSE, cexRow= 1.2,
                                           RowSideColors = ifelse( goont[ match( rownames(MMpval) ,names(goont))]=="BP" ,
                                                                   "#D55E00", ifelse( goont[ match( rownames(MMpval) ,names(goont))]=="CC" ,
                                                                                      "#009E73", "#0072B2")),
                                           ColSideColors = map2color(log10(sapply( GOlist , length)),
                                                                     pal = gray.colors(9, rev = T)))
          # dev.copy(which=a)
          # dev.off()
          # dev.off()
          print(dim(toPlot))
          
        }
        
        datTOplot <- Reduce( rbind , lapply( iid , function(ii){
          print(ii)
          datt <- datGO[[ii]]
          # select results for a given GO term, even if it's not a primary term in a given comparison
          datt <- datt$results[ match( GOtoPlot , datt$results$GO.ID )  , ]
          datt$dataset <- names(datGO)[ii]
          
          print( dim(datt) )
          
          ### add LFC or use log2Enrichment of enrichment
          if( !exists("datLFC") ) {
            datt$LFC <- datt[GOtoPlot , 'log2Enrichment']
          } else datt$LFC <- datLFC[GOtoPlot, ii]
          
          
          # make sore no terms are missing 
          datt$Ontology <- Ontology(GOTERM)[GOtoPlot]
          datt$Term <- Term(GOTERM)[GOtoPlot]
          datt$GO.ID <- GOtoPlot
          datt$Fisher[ is.na(datt$Fisher)] <-  1
          if( isTRUE(typeSig) ) datt$Ontology <- names(GOtoPlot)
          return(datt)
        }))
        
        
        # reorder levels
        if( exists("datOrder")) datTOplot$dataset <- factor( datTOplot$dataset,
                                     levels = datOrder)
        
        # reorder  GO IDs using the clustering order
        datTOplot$GO.ID <- factor( as.factor(datTOplot$GO.ID) ,
                                     levels= rownames(MMpval)[(gsea.order$rowInd)])
         

        # truncate labels
        if( isTRUE(labelTRUNK)) datTOplot$Term <- sapply( datTOplot$Term , str_trunc, width = 40)
        print(dim(datTOplot))
        
        # reorder Term using order of GO IDs from the clustering
        datTOplot$Term <- as.factor( datTOplot$Term)
        datTOplot$Term <-  factor( datTOplot$Term,
                                   levels=unique(datTOplot$Term[order(datTOplot$GO.ID)]), ordered=TRUE)
        
        # limit x axis
        datTOplot$log10.Fisher <- -log10(datTOplot$Fisher)
        datTOplot$log10.Fisher[ datTOplot$log10.Fisher>15] <- 15 
        # plot
        gg<- ggplot( datTOplot, aes(x=log10.Fisher, 
                                    y=Term ,
                                    fill=LFC ) ) + 
          scale_fill_continuous_divergingx(palette = "Spectral",
                                           na.value="grey",
                                           rev =T, mid = 0 ) +  
          scale_size(range = c(0, 15))+ ylab("GO terms") +
          # geom_vline(xintercept=0, linetype=2,lwd=2)+ # add vertical line for lfc=0
          geom_bar(stat = "identity",colour="black") + 
          facet_grid(cols = vars(dataset), 
                     # rows = vars(Ontology), 
                     scales = "free_y", space = "free")+
          theme_bw() + theme( text = element_text(size = 20 ),
                              axis.text.x = element_text(size = 16))
        
        return(gg)
      }
      
    # barplot for SPIA results     
      SPIA_barplot_clustered <-  function(iid=iids, 
                                          SPIAtoPlot= list(KEGG.spia_common,REACT.spia_common) ,
                                          SPIA_strictBkg_KEGG=SPIA_strictBkg_KEGG,
                                          SPIA_strictBkg_REACT= SPIA_strictBkg_REACT,
                                          labelTRUNK=T,
                                          typeSig=T)
      {
        ## * iids - a numeric vector, defining which elements of datGO list to plot as columns of the barplot
        ## * SPIAtoPlot - a list of 2 charchter vectors of 1) KEGG and 2) SPIA 
        ## * datSPIA - a dataframe with results of runSPIA()
        ## * oorder - a charachter vector, defining order in which datasets (in columns) appear on the plot
        ## * labelTRUNK - logical value, if labels should be shorten
        ## * typeSig - logical value, whether to split rows by Ontology type of universality of terms
        require(colorspace)
        library(GeneOverlap)
        
        ## combine reactome and kegg
        KFO.all_spia.both <- Reduce( rbind, c( SPIA_strictBkg_KEGG[iids] , 
                                               SPIA_strictBkg_REACT[iids] )  )
        
        # add database type 
        KFO.all_spia.both$dbase <- c( rep("KEGG", sum(sapply(SPIA_strictBkg_KEGG[iids], nrow))), 
                                      rep("REACTOME", sum(sapply(SPIA_strictBkg_REACT[iids], nrow))))
        KFO.all_spia.both$dataset <-c( rep.int( names(SPIA_strictBkg_KEGG)[iids] , 
                                                times= sapply(SPIA_strictBkg_KEGG[iids] , nrow)),
                                       rep.int( names(SPIA_strictBkg_REACT)[iids] , 
                                               times= sapply(SPIA_strictBkg_REACT[iids] , nrow)))
        # order columns
        KFO.all_spia.both$dataset <- factor( KFO.all_spia.both$dataset, 
                                             levels = names(SPIA_strictBkg_REACT)[iids])
        
        
        ### cluster for similarity of paths
        {
          # table of similarity
          Plist <- c( kegg_pathsList[ names(kegg_pathsList)%in% SPIAtoPlot[[1]]], 
                      reactPath_gName[ names(reactPath_gName)%in% SPIAtoPlot[[2]]])
          MM<- newGOM( Plist , Plist ,genome.size=10000)
          MMpval <- getMatrix(MM, name="pval")
          MMpval[MMpval==0] <- NA
          diag(MMpval) <- NA
          MMpval <- -log10(MMpval)
          
          ## plot
          rownames(MMpval)<- colnames(MMpval) <- make.unique( rownames(MMpval))
          toPlot <- MMpval  
          
          # pdf( width =  12 , height = 9, file="gseaKEGG.REACT.top3_pvalHeat.clust.pdf")
          # a<-dev.cur()
          # png( width =  1200 , height = 800 , file="gseaKEGG.REACT.top3_pvalHeat.clust.png")
          # dev.control("enable")
          gsea.order <- gplots::heatmap.2( toPlot,
                                           col = RColorBrewer::brewer.pal( name = "Greens", n=9), 
                                           trace = "none", symkey=T,
                                           na.color = "#00441B",
                                           dendrogram = "row",margins = c(2,20) ,
                                           labCol = FALSE, cexRow= 1.2,
                                           RowSideColors = ifelse( rownames(MMpval) %in% SPIAtoPlot[[1]] ,
                                                                   "royalblue" ,"violet" ),
                                           ColSideColors = map2color(log10(sapply( Plist , length)),
                                                                     pal = gray.colors(9, rev = T)))
          # dev.copy(which=a)
          # dev.off()
          # dev.off()
          print(dim(toPlot))
          
        }
        
        datTOplot <- KFO.all_spia.both[
          (KFO.all_spia.both$Name %in% SPIAtoPlot[[1]] & 
             KFO.all_spia.both$Name %in% names(kegg_pathsList)&
             KFO.all_spia.both$dbase=="KEGG") |
          (KFO.all_spia.both$Name %in% SPIAtoPlot[[2]] &
             KFO.all_spia.both$Name %in% names(reactPath_gName)&
             KFO.all_spia.both$dbase=="REACTOME"),]

        # # reorder levels
        # datTOplot$dataset <- factor( datTOplot$dataset,
        #                              levels = datOrder)
        datTOplot$Name <- factor( as.factor(datTOplot$Name) ,
                                  levels= rownames(MMpval)[(gsea.order$rowInd)])
        # truncate labels
        if( isTRUE(labelTRUNK)) levels(datTOplot$Name) <- sapply( levels(datTOplot$Name) , str_trunc, width = 40)
        print(dim(datTOplot))
        
        # # reorder Term using order of GO IDs from the clustering
        # datTOplot$Name <- as.factor( datTOplot$Name)
        # datTOplot$Name <-  factor( datTOplot$Name,
        #                            levels=unique(datTOplot$Name[order(datTOplot$GO.ID)]), ordered=TRUE)
        # 
        # limit x axis
        datTOplot$log10.pG <- -log10(datTOplot$pG)
        datTOplot$log10.pG[ datTOplot$log10.pG>10] <- 10 
        
        # plot
        gg<- ggplot(datTOplot, aes(
          x= log10.pG , 
          y= Name , 
          fill=Status)) + 
          scale_fill_manual( values = c("#E69F00","#56B4E9")) + 
          # scale_fill_colorblind() + 
          xlab("-log10(pG)")+
          geom_bar(stat = "identity") + ylab("pathways")+
          facet_grid(cols = vars(dataset), 
                     # rows = vars(dbase), 
                     scales = "free_y", space = "free")+
          theme_bw() + theme( text = element_text(size = 24, face = "bold" ),
                              axis.text.y = element_text(face = "plain" ),
                              axis.text.x = element_text(size = 16)) 
        
        return(gg)
      }
      
####======== 2D plots ==========

#### 2D plot for Robert's GO results:
## compares mean LFC (ctrl vs exprmnt) of GO terms in 2 conditions, 
## a term should be significant for at least one of the 2 conditions
  GO_2Dplot <- function( GOrb_list  , LFCdat , 
                          ids , dataName , contrasts , cnfdlvl=0.99 ){
     
     ## * GOrb_list - list of results for GO annotation, generated by sf.clusterGoByGeneset()
     ## * LFCdat - a numeric matrix or dataframe, containing mean LFC of GO gene members, 
     ## column names should equal to GOrb_list names
     ## * ids - 2 numeric values, ids of 2 datasets to be compared, will be used
     ## to select ids elements from GOrb_list and ids columns from LFCdat
     ## * dataName - character vector, descriptive names of 2 DEtests
     ## * contrasts - character vector, labels that differntiate 2 DEtests
     # x axis - first dataset, y-axis - second
     ## * cnfdlvl - a numeric value between 0 and 1 for the confidence interval
     # of a regression that fits mean LFCs of GO terms of 2 FSGS comparisons:
     # only terms outside of the conf.int will be labeled
     
     require( ggrepel )
     
     # get GO terms significant in at least one stage
     X <- union(GOrb_list[[ids[1]]]$results$GO.ID[which(GOrb_list[[ids[1]]]$results$Fisher< thrsh)],
                GOrb_list[[ids[2]]]$results$GO.ID[which(GOrb_list[[ids[2]]]$results$Fisher< thrsh )] )
     # # make sure ID is tested in both
     # X <- X[X %in% intersect(GOrb_list[[1]]$results$GO.ID, GOrb_list[[2]]$results$GO.ID)]
     # get mean LFC for the terms of interest
     GOrobert_LFC <- as.data.frame( LFCdat [ rownames( LFCdat ) %in% X , ids] )
     GOrobert_LFC[is.na(GOrobert_LFC)]<- 0 
     GOrobert_LFC <- GOrobert_LFC[rowSums(GOrobert_LFC!=0)>0,]
     # add term size
     # GOrobert_LFC$size <- sapply( 1:nrow(GOrobert_LFC), function(ii){
     #   length( All_GO$ensembl_gene_id[All_GO$go_id == rownames(GOrobert_LFC)[ii]])
     # } )
     GOrobert_LFC$FC.FC <- 10^(GOrobert_LFC[,2] - GOrobert_LFC[,1])
     GOrobert_LFC$pval.FC <- -1*log10(GOrb_list[[ids[2]]]$results$Fisher[match(  rownames(GOrobert_LFC) , GOrb_list[[ids[2]]]$results$GO.ID )] ) - 
       -1*log10(GOrb_list[[ids[1]]]$results$Fisher[ match( rownames(GOrobert_LFC), GOrb_list[[ids[1]]]$results$GO.ID )] )
     GOrobert_LFC$pval.4w <- GOrb_list[[ids[1]]]$results$Fisher[match(  rownames(GOrobert_LFC) , GOrb_list[[ids[1]]]$results$GO.ID )]
     GOrobert_LFC$pval.12w <- GOrb_list[[ids[2]]]$results$Fisher[match(  rownames(GOrobert_LFC) , GOrb_list[[ids[2]]]$results$GO.ID )]
     GOrobert_LFC$pval.4w[is.na(GOrobert_LFC$pval.4w)] <- 1
     GOrobert_LFC$pval.12w[is.na(GOrobert_LFC$pval.12w)] <- 1
     
     # # shorten the names
     # GOterm <- as.list(GOTERM)
     # GOrobert_LFC$Primary <- sapply( GOterm[rownames(GOrobert_LFC)], 
     #                                 function(X) X@Term )
     GOrobert_LFC$Primary <- goterms[ match( rownames(GOrobert_LFC), names(goterms))]
     
     # Create prediction interval data frame with upper and lower lines corresponding to sequence covering minimum and maximum of x values in original dataset
     colnames(GOrobert_LFC)[1:2 ]<- c("early","late")
     newx <- lm(  formula = GOrobert_LFC$late ~ GOrobert_LFC$early, )
     pred_interval <- predict(newx, interval="prediction", level = cnfdlvl)
     pred_interval <- as.data.frame(pred_interval)
     GOrobert_LFC <- cbind(GOrobert_LFC , pred_interval)
     pred_interval <- cbind(pred_interval , GOrobert_LFC[,1:6] )
     
     # 
     GOrobert_LFC$GOlabel <- ifelse( 
       GOrobert_LFC$late < GOrobert_LFC$lwr  |
         GOrobert_LFC$late > GOrobert_LFC$upr ,
       as.character(GOrobert_LFC$Primary), "")
     
     # add size
     GOrobert_LFC$size <- sapply( rownames(GOrobert_LFC), 
                                  function(X) length( unique( go2gName[[X]] )))
     # add type of GO
     GOrobert_LFC$GOtype <-  ifelse( 
       GOrobert_LFC$late < GOrobert_LFC$lwr  |
         GOrobert_LFC$late > GOrobert_LFC$upr ,
       as.character(  Ontology(GOTERM)[rownames(GOrobert_LFC)]), " ")

# plot
    gg <-  ggplot( GOrobert_LFC , aes(x= early ,y= late) ) +  
       geom_point(aes( size = size ,  color= GOtype ), shape=21, stroke = 2) + 
       scale_color_manual( values = c("BP" ="#D55E00", "CC" ="#009E73", 
                                      "MF" ="#0072B2", " " = "grey" ) )+
       geom_ribbon( data=pred_interval, aes(ymin = lwr, ymax = upr), fill = "blue", alpha = 0.1) + 
       # geom_line( data=GOrobert_LFC, aes(x = fit, y = upr), colour = "darkgrey")+
       # geom_line( data=GOrobert_LFC, aes(x = fit, y = lwr), colour = "darkgrey")+
       geom_hline(yintercept=0, linetype="dashed", color = "black")+
       geom_vline(xintercept=0, linetype="dashed", color = "black")+
       scale_size(range = c(1, 30)) +  # adjust range of bubbles sizes
       geom_text_repel( aes(label = GOlabel ,colour=GOtype), 
                        max.overlaps=75, cex=6 ,
                        max.time	=100) + #ensure that labels are not jammed together
       xlab(paste("mean log2FC of GO at", contrasts[1],"weeks", sep = " ")) + 
       ylab(paste("mean log2FC of GO at", contrasts[2],"weeks", sep = " ")) + 
       ggtitle(paste( "activity changes in FSGS-related GO terms\n",
                      dataName, sep = "")) + theme_bw()+
       # coord_cartesian( xlim = c(-0.2, max(GOrobert_LFC$Wt1het.del_4w)) ,
       #                  ylim = c(-0.2, max(GOrobert_LFC$Wt1het.del_12w)) ) +
       theme( text = element_text(size = 24))
    return(gg)
   }   

#### 2-way ^plot of Robert's GO results: 
## mean LFC VS Pvalue
      
    GO_LFC.PVALplot <- function( GOrbrt=RobertGO.50 , LFCvec = LFCvec,
                                 dataName ,trunc =50 ){
      
      ## * GOrbrt - results of GO annotation, generated by sf.clusterGoByGeneset()
      ## * LFCdat - a numeric vector, containing mean LFC of GO gene members
      
      require( ggrepel )
      require(stringr)
      # get mean LFC for the terms of interest
      LFCvec <-  LFCvec
      # 
      GOrobert_LFC <- GOrbrt$results
      # GOrobert_LFC <- GOrobert_LFC[ GOrobert_LFC$Is.primary=="TRUE", ]
      GOrobert_LFC$LFC <- LFCvec[ match( GOrobert_LFC$GO.ID, names(LFCvec))]
      # 
      GOrobert_LFC$GOlabel <- ifelse( GOrobert_LFC$Primary.Fisher.adj < thrsh & 
                                        GOrobert_LFC$Is.primary=="TRUE" &
                                        abs(GOrobert_LFC$LFC)>=0.01 ,
                                      as.character(GOrobert_LFC$Term), "")
      
      # add size
      GOrobert_LFC$size <- sapply( rownames(GOrobert_LFC), 
                                   function(X) length( unique( go2gName[[X]] )))
      # # add type of GO
      # GOrobert_LFC$GOtype <- as.character(  Ontology(GOTERM)[rownames(GOrobert_LFC)])
      # 
      # transform pval
      GOrobert_LFC$Fisher.m.log10 <- -log10(GOrobert_LFC$Fisher)
      GOrobert_LFC$Padj.m.log10 <- -log10(GOrobert_LFC$Primary.Fisher.adj)
      # limit span of pval
      GOrobert_LFC$Fisher.m.log10[GOrobert_LFC$Fisher.m.log10>7] <- 7
      GOrobert_LFC$Padj.m.log10[GOrobert_LFC$Padj.m.log10>7] <- 7
      GOrobert_LFC$LFC[GOrobert_LFC$LFC>0.05] <- 0.05
      GOrobert_LFC$LFC[GOrobert_LFC$LFC<-0.05] <- -0.05
      
      
      # trunc
      if( trunc>0) GOrobert_LFC$GOlabel <- sapply(GOrobert_LFC$GOlabel , str_trunc, trunc)
     
       # plot
      gg <- ggplot( GOrobert_LFC , aes(x= LFC ,y= Fisher.m.log10 ) ) +  
        geom_point(aes( size= Is.primary, fill=LFC), shape=21, stroke = 0) + 
        scale_fill_distiller( palette = "PiYG" , 
                              limit=c(-0.03, 0.03), na.value = "#9E0142")+
        geom_hline(yintercept=-log10(thrsh), linetype="dashed", color = "black")+
        geom_vline(xintercept=-0.01, linetype="dashed", color = "black")+
        geom_vline(xintercept=0.01, linetype="dashed", color = "black")+
        geom_text_repel( aes(label = GOlabel ,colour=Ontology), max.overlaps=100, cex=6 ,
                         max.time	=50) + #ensure that labels are not jammed together
        scale_color_manual( values = c("BP" ="#D55E00", "CC" ="#009E73", 
                                       "MF" ="#0072B2" ) )+
        
        xlab( "mean LFC of GO") + ylab( "-log10( Fisher.pval )") + 
        ggtitle(paste( "2-way plot of GO terms\n",
                       dataName, sep = "")) + theme_bw()+
        # coord_cartesian( xlim = c(-0.2, max(GOrobert_LFC$Wt1het.del_4w)) ,
        #                  ylim = c(-0.2, max(GOrobert_LFC$Wt1het.del_12w)) ) +
        theme( text = element_text(size = 24))
      return(gg)
    }   
      


# ####======== Circular chord plot ==========
# # description on http://wencke.github.io/ 
# 
# # # choose pathways to plot
# # path_sel <- read.table(sep = "\t","/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/SPIA/SPIAKEGGREACTOMEOVERLAP_2018.csv")
# # path_selKEGG <- as.character(path_sel[which(path_sel[,2]=="KEGG"),1])
# # path_selREACT <- as.character(path_sel[which(path_sel[,2]=="REACTOME"),1])
# 
# 
# ### function to make a chord plot given results of SPIA and DEseq2
# chord_plotLFC <- function( expr1=WT1_4w , expr2=WT1_12w , sortBy = 2, 
#                            DEseq2qvalthresh=0.01 , path_sel = paths )
#   {
#   ### paths - names of kegg and/or reactome pathways to plot 
#   
#   
#   # require(graphite)
#   require("qdapTools")
#   
#   # load necessary functions 
#   source("/home/tim_nevelsk/PROJECTS/myCode/Helper_chordPlot.R")
#   
#   ### get pathway gene sets
#   {  # parepare Reactome
#     {
#       # read the file
#       reactPath <- fgsea::gmtPathways( "/media/tim_nevelsk/WD_tim/ANNOTATIONS/Pathways/Reactome/ReactomePathways.gmt")
#       # convert human to mouse IDS
#       
#       fun_homoTOmouse <- function(gns){
#         egs <- mapIds(org.Hs.eg.db, gns, "ENTREZID","SYMBOL")
#         mapped <- select(Orthology.eg.db, egs, "Mus.musculus","Homo.sapiens")
#         mapped$MUS <- mapIds(org.Mm.eg.db, as.character(mapped$Mus.musculus), "SYMBOL", "ENTREZID")
#         mapped$MUS.ens <- mapIds(org.Mm.eg.db, as.character(mapped$Mus.musculus), "ENSEMBL", "ENTREZID")
#         return(mapped)
#       }
#       
#       options(connectionObserver = NULL)
#       library( "org.Hs.eg.db" )
#       library("org.Mm.eg.db")
#       library(Orthology.eg.db)
#       
#       humanGen <- unique (Reduce( c , reactPath) )
#       mouseGen <- fun_homoTOmouse(humanGen)
#       mouseGen$HOMO <- rownames(mouseGen)
#       
#       reactPath_gENSID <- lapply( reactPath , function(x){
#         X <- as.character( na.omit( unique(mouseGen$MUS.ens[ match( x, mouseGen$HOMO)])))
#       } )
#       reactPath_gName <- lapply( reactPath , function(x){
#         X <- as.character( na.omit( unique(mouseGen$MUS[ match( x, mouseGen$HOMO)])))
#       } )
#       
#       # remove untranslated names
#       reactPath_gName <- lapply(reactPath_gName, function(X) X <- X[X!="NULL"])
#       reactPath_gENSID <- lapply(reactPath_gENSID, function(X) X <- X[X!="NULL"])
#       
#     }
#     
#     # parepare KEGG
#     {
#       ## download directly from KEGG: 341 paths
#       library(EnrichmentBrowser)
#       options(connectionObserver = NULL)
#       
#       library(KEGGgraph)
#       keggPath_gName <- EnrichmentBrowser::getGenesets( "mmu", db ="kegg",
#                                                         gene.id.type = "SYMBOL" )
#       
#       # names(keggPath_gName) <- gsub( "_"," " ,gsub("^[^_]*_" , "",  names(keggPath_gName)))
#       
#     }}
# 
#   
#   ## load graphite pathway data
#   # load("/data/user/tpadvits/PROJECTS/MARTIN_PJ/RNAseq_analysis/NPHS2/SPIA/SPIApathInfo.02.05.19.rda")
#   
#     SEL_LFC <- c(keggPath_gName[ gsub( "_"," " ,
#                                        gsub("^[^_]*_" , "",  names(keggPath_gName))) %in% paths ] ,
#                  reactPath_gName[ names(reactPath_gName) %in% paths ]) 
#     
#   
#   ### LFC of genes from Selected KEGG and Reactome pathways 
#   # convert lists to df
#     SEL_LFC_df <- list2df(SEL_LFC)
# 
# 
#   ### add LFC
#   if( "lfc_shrunk" %in% colnames(expr1) ) { 
#     lfc_colName <- "lfc_shrunk"  } else  lfc_colName <- "log2FoldChange"
#   
#     SEL_LFC_df$w4lfc <- expr1[[lfc_colName]][ match(SEL_LFC_df$X1,expr1$gName)]
#     SEL_LFC_df$w4qval <-  expr1$padj[ match(SEL_LFC_df$X1,expr1$gName) ]
#     SEL_LFC_df$w12lfc <- expr2[[lfc_colName]][ match(SEL_LFC_df$X1,expr2$gName)]
#     SEL_LFC_df$w12qval <-  expr2$padj[ match(SEL_LFC_df$X1,expr2$gName)]
#     SEL_LFC_dfSig <- subset( SEL_LFC_df, SEL_LFC_df$w4qval< DEseq2qvalthresh |
#                                SEL_LFC_df$w12qval< DEseq2qvalthresh )
#     # SEL_LFC_dfSig$X2 <- paste( SEL_LFC_df$X2,"REACTOME", sep = "^")
#     
#    
#  
#   
#   colnames(SEL_LFC_dfSig) <- c("gName","pathway", "lfc1" ,"qval1" , "lfc2" ,"qval2" )
#   X_4w <- SEL_LFC_dfSig[,c(2,1,3)]
#   X_12w <- SEL_LFC_dfSig[,c(2,1,5)]
#     
#   colnames(X_4w) <- colnames(X_12w) <-  c("term","genes","logFC")
#   chord_4w <- chord_datX(X_4w)
#   chord_12w <- chord_datX(X_12w)
#   
#   
#   #* Add an warning if sortBy not 1 or 2
#   if ( !(sortBy%in%c(1,2)) ){
#       warning( "'sortBy' must be 1 or 2. sortBy is set to 1" )
#     sortBy <- 1
#     }
#     
#   if ( sortBy ==1) {
#     chord_4w <- chord_4w[order(-chord_4w[,ncol(chord_4w)]),]
#     chord_12w <- chord_12w[rownames(chord_4w),]
#   } else if (sortBy ==2) {
#     chord_12w <- chord_12w[order(-chord_12w[,ncol(chord_12w)]),]
#     chord_4w <- chord_4w[rownames(chord_12w),]
#   }
# 
#   chord_4w12w <- cbind.data.frame(chord_4w,logFC=chord_12w[,ncol(chord_12w)])
#   
#   pp<- GOChordX( chord_4w12w, space = 0.01, gene.order = 'none',process.label=8	,
#             gene.space = 0.25, gene.size = 4 , nlfc=2 , border.size =0.02,  
#             lfc.min=min(c(chord_4w[,ncol(chord_4w)],chord_12w[,ncol(chord_12w)])), 
#             lfc.max=max(c(chord_4w[,ncol(chord_4w)],chord_12w[,ncol(chord_12w)])) )
#   print(pp)
#   
#   return( SEL_LFC_dfSig)
# }
# 
# 
#   
# 
# 
# ####======== Heatmaps of important genes ==========
# # gene targets of WT1hetdel mouse
# WT1_closGene <- read.table( sep = "\t" , header = F , file="/data/user/tpadvits/PROJECTS/MARTIN_PJ/DATA_integ/WT1_Tead1/coocurance/WT1_nearestg.bed")
# WT1_closGene_sig <- unique (sub("\\..*","",WT1_closGene$V10[which(abs(WT1_closGene$V13)<=2000)]))
# WT1_closGene_sig_genes <- unique(tx2gene$external_gene_name[which(tx2gene$ensembl_transcript_id%in%(WT1_closGene_sig))])
# 
# 
# ### interesting SPIA pathways
# spia_results <- list.files(pattern = "w_SPIA_","/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/SPIA/",full.names = T)
# ephr_paths_react <- do.call( union , Filter(length,lapply(spia_results, function(x) grep("ephr",readRDS(x)$Name,value = T,ignore.case = T))))
# coll_paths_react <- do.call( union , Filter(length,lapply(spia_results, function(x) grep("collag",readRDS(x)$Name,value = T,ignore.case = T))))
# 
# 
# ###============= heatmap of the expression across replicates
# # read expression data
# expression_dat <- read.table(header = T, stringsAsFactors = F, row.names = 1, "/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/DE_analysis/expr_IntExon_genes_tpmSized.csv")
# colnames(expression_dat) =c("w4_het.del","w4_wt","w4_wt","w4_het.del","w4_wt","w4_het.del","w12_wt","w12_wt","w12_het.del","w12_het.del","w12_het.del","w12_wt")
# expression_dat$uniprotID <- tx2gene$uniprot_gn[match(rownames(expression_dat),tx2gene$ensembl_gene_id)]
# expression_dat$gName <- tx2gene$external_gene_name[match(rownames(expression_dat),tx2gene$ensembl_gene_id)]
# expression_dat_filt <- expression_dat[which(expression_dat$gName %in% WT1_4w$gName),]
# 
# ### Ephrin and collagene pathways
# {   
#   ephr_react_TPM <- lapply( ephr_paths_react, function(x) { pathGenes_entrID <- sub("UNIPROT:","",nodes(path_REACT[[x]]))
#     l2FC <- expression_dat_filt[which(expression_dat_filt$uniprotID%in%pathGenes_entrID),1:12]
#     rownames(l2FC) <-  expression_dat_filt$gName[which(expression_dat_filt$uniprotID%in%pathGenes_entrID)]
#     return(l2FC) })
#   ephr_react_TPM <- do.call(rbind, ephr_react_TPM)
#   ephr_react_TPM <- ephr_react_TPM[!duplicated(ephr_react_TPM),]
#   
#   collag_react_TPM <- lapply( coll_paths_react, function(x) { pathGenes_entrID <- sub("UNIPROT:","",nodes(path_REACT[[x]]))
#     l2FC <-  expression_dat_filt[which(expression_dat_filt$uniprotID%in%pathGenes_entrID),1:12]
#     rownames(l2FC) <-  expression_dat_filt$gName[which(expression_dat_filt$uniprotID%in%pathGenes_entrID)]
#     return(l2FC) })
#   collag_react_TPM <- do.call(rbind, collag_react_TPM)
#   collag_react_TPM <- collag_react_TPM[!duplicated(collag_react_TPM),]
#   
# } 
# ### collagene genes 
# { collagene_TPM <- read.table(sep = "\t", stringsAsFactors = F, row.names = 1, header = T, "/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/expr_heatmaps/collagens_expr/collagen_tpm.csv")
#   collagene_TPM <- collagene_TPM[which(rownames(collagene_TPM)%in%expression_dat_filt$gName),]
#   collagene_TPM_ave <- cbind(rowMeans(collagene_TPM[,1:3]),rowMeans(collagene_TPM[,4:6]),rowMeans(collagene_TPM[,7:9]),rowMeans(collagene_TPM[,10:12]))
#   colnames(collagene_TPM_ave) <- c("w4_wt","w4_het.del","w12_wt","w12_het.del")
# }
# 
# library(ggplot2)
# library(reshape)
# #XX <- t(scale(t(ephr_react_LFC),center = F))
# XX <- collagene_TPM_ave[which(rowMeans(collagene_TPM_ave) > 1),]
# XX <- t(scale(t(XX[,1:4]),center = F))
# ord <- hclust( dist(XX, method = "euclidean"), method = "ward.D" )$order
# XX.m <- melt( as.matrix( XX ) )
# XX.m$X1 <- factor( XX.m$X1, levels = rownames(XX)[ord], labels =  rownames(XX)[ord] )
# 
# 
# ggplot(XX.m, aes(X2, X1)) + 
#   geom_tile(aes(fill = value), colour = "darkgrey") + theme_bw() + theme( panel.grid.major = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
#                                                                          panel.border = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
#                                                                                              axis.text.x = element_text(face="bold", angle=45,vjust=0.5,size=12),axis.text.y = element_text(size=10,face="italic",family = "arial"))+ 
#   # geom_text(aes(label = round(padjust, 4))) +    axis.text.x = element_text(face="bold", angle=0,vjust=0.5,size=11),axis.text.y = element_text(size=11,face="italic",family = "arial"))
#   scale_fill_gradient2( midpoint = 0.8, low="white" ,mid = "yellow", high = "red", na.value = "grey",name="scaled tpm of gene expression") 
#   # scale_fill_manual(values=my_palette, breaks=c(levels(as.factor(XX$value))[seq(1, 160, by=10)],2.979),name="log2 Fold.Change",na.value = "red")
# 
# ###================ LFC heatmap 
# ### Ephrin and collagene pathways
# { 
#   ephr_react_LFC <- lapply( ephr_paths_react, function(x) { pathGenes_entrID <- sub("UNIPROT:","",nodes(path_REACT[[x]]))
#   l2FC <-  cbind( WT1_4w$log2FoldChange[which(WT1_4w$uniprotID%in%pathGenes_entrID)] , WT1_12w$log2FoldChange[which(WT1_12w$uniprotID%in%pathGenes_entrID)] )
#   rownames(l2FC) <-  WT1_4w$gName[which(WT1_4w$uniprotID%in%pathGenes_entrID)]
#   return(l2FC) })
#   ephr_react_LFC <- do.call(rbind, ephr_react_LFC)
#   ephr_react_LFC <- ephr_react_LFC[!duplicated(ephr_react_LFC),]
#   
#   collag_react_LFC <- lapply( coll_paths_react, function(x) { pathGenes_entrID <- sub("UNIPROT:","",nodes(path_REACT[[x]]))
#   l2FC <-  cbind( WT1_4w$log2FoldChange[which(WT1_4w$uniprotID%in%pathGenes_entrID)] , WT1_12w$log2FoldChange[which(WT1_12w$uniprotID%in%pathGenes_entrID)] )
#   rownames(l2FC) <-  WT1_4w$gName[which(WT1_4w$uniprotID%in%pathGenes_entrID)]
#   return(l2FC) })
#   
#   collag_react_LFC <- do.call(rbind, collag_react_LFC)
#   collag_react_LFC <- collag_react_LFC[!duplicated(collag_react_LFC),]
#   
#   colnames(ephr_react_LFC)   <- colnames(collag_react_LFC) <- c("week 4", "week 12")
#   }
# 
# ### collagene genes 
# collagene_l2FC <-  cbind( WT1_4w[which(WT1_4w$gName%in%rownames(collagene_TPM)),c("log2FoldChange","gName")] , WT1_12w[which(WT1_12w$gName%in%rownames(collagene_TPM)),c("log2FoldChange","gName")] )
# rownames(collagene_l2FC) <- collagene_l2FC$gName
# collagene_l2FC <- collagene_l2FC[,c(1,3)]
# colnames(collagene_l2FC) <- c("week4","week12")
# 
# 
# library(ggplot2)
# library(reshape)
# #XX <- t(scale(t(ephr_react_LFC),center = F))
# XX <- collagene_l2FC
# #XX <- t(scale(t(collagene_l2FC),center = F))
# ord <- hclust( dist(XX, method = "euclidean"), method = "ward.D" )$order
# XX.m <- melt( as.matrix( XX ) )
# XX.m$X1 <- factor( XX.m$X1, levels = rownames(XX)[ord], labels =  rownames(XX)[ord] )
# 
# 
# ggplot(XX.m, aes(X2, X1)) + 
#   geom_tile(aes(fill = value), colour = "darkgrey") + theme_bw() + theme( panel.grid.major = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
#                                                                           panel.border = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
#                                                                           axis.text.x = element_text(face="bold", angle=45,vjust=0.5,size=12),axis.text.y = element_text(size=10,face="italic",family = "arial"))+ 
#   # geom_text(aes(label = round(padjust, 4))) +    axis.text.x = element_text(face="bold", angle=0,vjust=0.5,size=11),axis.text.y = element_text(size=11,face="italic",family = "arial"))
#   scale_fill_gradient2( low="blue" ,mid = "white", high = "red", na.value = "grey",name="scaled log2 fold change\nof gene expression") 
# # scale_fill_manual(values=my_palette, breaks=c(levels(as.factor(XX$value))[seq(1, 160, by=10)],2.979),name="log2 Fold.Change",na.value = "red")
# 
#  