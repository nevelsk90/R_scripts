#======== Circular chord plot ==========
# description on http://wencke.github.io/ 

# # choose pathways to plot
# path_sel <- read.table(sep = "\t","/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/SPIA/SPIAKEGGREACTOMEOVERLAP_2018.csv")
# path_selKEGG <- as.character(path_sel[which(path_sel[,2]=="KEGG"),1])
# path_selREACT <- as.character(path_sel[which(path_sel[,2]=="REACTOME"),1])


### function to make a chord plot given results of SPIA and DEseq2

chord_plotLFC <- function( expr1 , expr2 ,  
                           path_selKEGG = path_selKEGG,  sortBy = 2,  
                           path_selREACT = path_selREACT, DEseq2qvalthresh=0.01,
                           sizeP , filePref )
{
  
  require(graphite)
  require(ggplot2)
  require("qdapTools")
  require( path_REACT )
  require( path_KEGG )
  
  # load necessary functions 
  source("/home/tim_nevelsk/PROJECTS/myCode/Helper_chordPlot.R")
  
  # load graphite pathway data
  load("/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_NPHS2/SPIApathInfo.25.05.21.Rdata")
  
  SELreact_LFC <- lapply( path_selREACT , function(x) tryCatch( nodes(path_REACT[[x]]),error=function(e) numeric() ) )
  names(SELreact_LFC) <- path_selREACT
  ### LFC of genes from Selected KEGG and Reactome pathways 
  # convert lists to df
  SELreact_LFC_df <- list2df(SELreact_LFC)

  if (nrow(SELreact_LFC_df) !=0 ) { 
  # convert specific IDs to gene names
  SELreact_LFC_df[,1] <- tx2gene$external_gene_name[match(sub("ENTREZID:","",SELreact_LFC_df[,1]),tx2gene$entrezgene_id)] 

  # remove duplicate rows and NA genes 
  SELreact_LFC_df <- SELreact_LFC_df[!duplicated(SELreact_LFC_df),]
  SELreact_LFC_df <- SELreact_LFC_df[!is.na(SELreact_LFC_df[,1]),]
  
  # add LFC
  SELreact_LFC_df$w4lfc <- expr1$log2FoldChange[ match(SELreact_LFC_df$X1,expr1$gName)]
  SELreact_LFC_df$w4qval <-  expr1$padj[ match(SELreact_LFC_df$X1,expr1$gName) ]
  SELreact_LFC_df$w12lfc <- expr2$log2FoldChange[ match(SELreact_LFC_df$X1,expr2$gName)]
  SELreact_LFC_df$w12qval <-  expr2$padj[ match(SELreact_LFC_df$X1,expr2$gName)]
  SELreact_LFC_dfSig <- subset( SELreact_LFC_df, SELreact_LFC_df$w4qval< DEseq2qvalthresh |  SELreact_LFC_df$w12qval< DEseq2qvalthresh )
  SELreact_LFC_dfSig$X2 <- paste( SELreact_LFC_dfSig$X2,"REACTOME", sep = "^")

  } else SELreact_LFC_dfSig <- NULL
  
 SELkegg_LFC <- lapply( path_selKEGG , function(x)  tryCatch(nodes(path_KEGG[[x]]),error=function(e) numeric() ) )
  names(SELkegg_LFC) <- path_selKEGG 
  ### LFC of genes from Selected KEGG and Reactome pathways 
  # convert lists to df
  SELkegg_LFC_df <-  list2df(SELkegg_LFC)
  
  if ( nrow(SELkegg_LFC_df) != 0 ) { 
  # convert specific IDs to gene names
  SELkegg_LFC_df[,1] <- tx2gene$external_gene_name[match(sub("ENTREZID:","",SELkegg_LFC_df[,1]),tx2gene$entrezgene_id)]
  
  # remove duplicate rows and NA genes 
  SELkegg_LFC_df <- SELkegg_LFC_df[!duplicated(SELkegg_LFC_df),]
  SELkegg_LFC_df <- SELkegg_LFC_df[!is.na(SELkegg_LFC_df[,1]),]
  
  # add LFC
  SELkegg_LFC_df$w4lfc <- expr1$log2FoldChange[match(SELkegg_LFC_df$X1,expr1$gName)]
  SELkegg_LFC_df$w4qval <-  expr1$padj[match(SELkegg_LFC_df$X1,expr1$gName)]
  SELkegg_LFC_df$w12lfc <- expr2$log2FoldChange[match(SELkegg_LFC_df$X1,expr2$gName)]
  SELkegg_LFC_df$w12qval <-  expr2$padj[match(SELkegg_LFC_df$X1,expr2$gName)]
  SELkegg_LFC_dfSig <- subset(SELkegg_LFC_df, SELkegg_LFC_df$w4qval< DEseq2qvalthresh |  SELkegg_LFC_df$w12qval< DEseq2qvalthresh )
  SELkegg_LFC_dfSig$X2 <- paste(SELkegg_LFC_dfSig$X2,"KEGG", sep = "^")
  } else SELkegg_LFC_dfSig <- NULL
  
  
  # # choose significant patwhays
  # all_SPIA_REACT_sig <- lapply(all_SPIA_REACT, function(x) subset(x, x$pGFdr< pGFDRthresh))
  # all_SPIA_KEGG_sig <- lapply(all_SPIA_KEGG, function(x) subset(x, x$pGFdr< pGFDRthresh))
  # 
  
  SELboth_LFC_dfSig <- rbind( SELreact_LFC_dfSig, SELkegg_LFC_dfSig )
  colnames(SELboth_LFC_dfSig) <- c("gName","pathway", "lfc1" ,"qval1" , "lfc2" ,"qval2" )
  X_4w <- SELboth_LFC_dfSig[,c(2,1,3)]
  X_12w <- SELboth_LFC_dfSig[,c(2,1,5)]
  
  colnames(X_4w) <- colnames(X_12w) <-  c("term","genes","logFC")
  chord_4w <- chord_datX(X_4w)
  chord_12w <- chord_datX(X_12w)
  
  
  #* Add an warning if sortBy not 1 or 2
  if ( !(sortBy%in%c(1,2)) ){
    warning( "'sortBy' must be 1 or 2. sortBy is set to 1" )
    sortBy <- 1
  }
  
  if ( sortBy ==1) {
    chord_4w <- chord_4w[order(-chord_4w[,ncol(chord_4w)]),]
    chord_12w <- chord_12w[rownames(chord_4w),]
  } else if (sortBy ==2) {
    chord_12w <- chord_12w[order(-chord_12w[,ncol(chord_12w)]),]
    chord_4w <- chord_4w[rownames(chord_12w),]
  }
  
  chord_4w12w <- cbind.data.frame(chord_4w,logFC=chord_12w[,ncol(chord_12w)])
  
  pp<- GOChordX( chord_4w12w, space = 0.01, gene.order = 'none',process.label=8	,
                 gene.space = 0.25, gene.size = 4 , nlfc=2 , border.size =0.02,  
                 lfc.min=min(c(chord_4w[,ncol(chord_4w)],chord_12w[,ncol(chord_12w)])), lfc.max=max(c(chord_4w[,ncol(chord_4w)],chord_12w[,ncol(chord_12w)])) )
  
  pdf( height = sizeP, width = sizeP, file= paste(filePref,".pdf",sep = ""))
  print(pp)
  dev.off()
  
  return( SELboth_LFC_dfSig)
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
