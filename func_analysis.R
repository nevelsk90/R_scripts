#=================================================================#
# ========== FUNCTIONAL ANNOTATION AND PATHWAY ANALYSIS ========= =
#=================================================================#
####======== GO ANALYSIS ============####
################# topGO analysis################ # 
library("topGO", lib.loc="/data/library/L-3.3.1")
# define gene set and gene universe
  {
    universe=rownames(txi$counts)
    geneset=rownames(DEseq_RUVsCon_l2FC_cond_clust[which(DEseq_RUVsCon_l2FC_cond_clust[,3]=="4"),])
    
    ### create named list that contain all gene members for each GO term
    mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "mmusculus_gene_ensembl",host = 'ensembl.org')
    All_GO=getBM(attributes=c("ensembl_transcript_id","ensembl_gene_id","go_id"), mart=mart)
    GO_list=unstack(All_GO[,c(2,3)]) 

    ### FUNCTION for topGO analysis 
    GoBytopGO_weight01 = function( GO_list, geneset, universe, min.genes=5, topGOalgorythm="weight01Count", testtype=GOFisherTest) 
      {
  require("topGO")
  
  ### create a gene list
  # first make a factor that is 1 if the probeset is "interesting" and 0 otherwise
  geneList <- factor(as.integer (All_ID %in% geneset))
  # name the factor with the probeset names
  names (geneList) <- All_ID
  
  # create topGOdata object for each category of GO
  topGO_BP=new("topGOdata",ontology='BP',description='topGOdata_BP', allGenes=geneList, annot = annFUN.GO2genes, GO2genes=GO_list, nodeSize = min.genes)
  topGO_MF=new("topGOdata",ontology='MF',description='topGOdata_MF', allGenes=geneList, annot = annFUN.GO2genes, GO2genes=GO_list, nodeSize = min.genes)
  topGO_CC=new("topGOdata",ontology='CC',description='topGOdata_CC', allGenes=geneList, annot = annFUN.GO2genes, GO2genes=GO_list, nodeSize = min.genes)
  
  #### test statistics definition ####
  #### OBS! when only a list of interesting genes is provided, the user can use only  
  #### tests statistics that are based on gene counts, like Fisher's exact test.
  test_stat <- new(topGOalgorythm, testStatistic = testtype, name = paste(topGOalgorythm,"_",testtype@generic[1]))
  
  ### test overrepresentation
  resultBP <- getSigGroups(topGO_BP, test_stat)
  resultMF <- getSigGroups(topGO_MF, test_stat)
  resultCC <- getSigGroups(topGO_CC, test_stat)
  
  #### create summary table from results with all GO terms ###
  # tables contain all GO terms used by each GO object, specified by length(usedGO(object = topGO_BP))
  gotable_BP <- GenTable(topGO_BP, weight01 = resultBF, ranksOf = "classic", topNodes = length(usedGO(object = topGO_BP)),numChar=100)
  gotable_MF <- GenTable(topGO_MF, weight01 = resultMF, ranksOf = "classic", topNodes = length(usedGO(object = topGO_MF)),numChar=100)
  gotable_CC <- GenTable(topGO_CC, weight01 = resultCC, ranksOf = "classic", topNodes = length(usedGO(object = topGO_CC)),numChar=100)
  
  # remove topGO with p-val=1 (don't contain any gene of interest) and add type of the GO terms 
  gotable_BP=cbind(type="BP",subset(gotable_BP,gotable_BP$weight01<1))
  gotable_MF=cbind(type="MF",subset(gotable_MF,gotable_MF$weight01<1))
  gotable_CC=cbind(type="CC",subset(gotable_CC,gotable_CC$weight01<1))
  
  # combine results for different GO categories in one table
  gotable=rbind(gotable_BP,gotable_MF,gotable_CC)
  gotable_pval0.01=subset(gotable,gotable$weight01<0.01)
  return(gotable)
}

# apply function
  DEseq_RUVsCon_gOnly_clust1.2_topGO_pval0.01=GoBytopGO_weight01(GO_list,geneset,universe)

## plot GO hierarchy
showSigOfNodes(so4ko12_GO, score(resultBF_so4ko12), firstSigNodes = 5, useInfo = "def")}

################# GO annotation with Robert function ################ # 
  {
  ### prepare gene set by annotating ensembleID with entrezID
    options(connectionObserver = NULL)
    
  source("/home/tim_nevelsk/PROJECTS/myCode/justGO_ROBERT.R")  
  library("org.Mm.eg.db")
    
  # use Robert's function to create sparse binary matrix indicating membership of genes (columns) in GO terms (rows)
  gomatrix=sf.createGoMatrix()
  
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
  tx2gene <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name","transcript_biotype","gene_biotype","entrezgene_id"), mart = mart)
  # tx2gene <- dplyr::rename(tx2gene, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id)
  
  # define background gene set
  universe=colnames(gomatrix)
  
  # prepare gene set, convert ensembleIDs to entrezIDs
  geneset=rownames()
  geneset=unique(tx2gene[which(tx2gene$ensembl_transcript_id%in%geneset),6])
  geneset=geneset[-1]
  geneset=as.character(geneset)
  
  # apply Robert's function that given a sparse matrix of GO terms (columns = genes, rows = GO terms)
  # a geneset of interest and a background set of genes (universe)
  # carry out clustering with members diverging by at most cut_max genes, and do enrichment testing.
  # Note, multiplicity adjustment is performed for the representative terms only.
  dbWT1_RobertGO=sf.clusterGoByGeneset(gomatrix,geneset,universe)
  dbWT1_RobertGO_qval0.1=subset(dbWT1_RobertGO$results, Enrichment>0 & Is.primary==TRUE & Primary.Fisher.adj<0.1)
  dbWT1_RobertGO_qval0.1=dbWT1_RobertGO_qval0.1[,4:13]
  
  #par(mar=par("mar") + c(3,18,3,3), xpd=F)
  margin_x=c(1,1,1,1)
  sf.plotResults(x)
}

#### make 2D plot of early vs late FSGS
  { # read data
  WT1_4w=read.table("/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/DE_analysis/DEexonIntron_4w_glvl.csv")
  WT1_12w=read.table("/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/DE_analysis/DEexonIntron_12w_glvl.csv")
  RobertGO_4w <- readRDS(file="/data/user/tpadvits/PROJECTS/MARTIN_PJ/Chipseq_analysis/WT1/Diff_bind/ExprBind_integration/RobertGO_4w.rds")
  RobertGO_12w <- readRDS(file="/data/user/tpadvits/PROJECTS/MARTIN_PJ/Chipseq_analysis/WT1/Diff_bind/ExprBind_integration/RobertGO_12w.rds")
  
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "mmusculus_gene_ensembl",host = "http://aug2017.archive.ensembl.org")
  geneID2entrez=biomaRt::getBM(attributes = c("ensembl_gene_id","entrezgene"), mart = mart)
  RobertGO_sigUni <- union(RobertGO_4w$GO.ID[which(RobertGO_4w$Primary.Fisher.adj<0.1)],RobertGO_12w$GO.ID[which(RobertGO_12w$Primary.Fisher.adj<0.1)])
  X <- union(RobertGO_4w$Term[which(RobertGO_4w$Primary.Fisher.adj<0.1)],RobertGO_12w$Term[which(RobertGO_12w$Primary.Fisher.adj<0.1)])
  
  library("plyr")
  source("/data/user/tpadvits/PROJECTS/MARTIN_PJ/DATA_integ/justGO_ROBERT.R")
  gomatrix <- sf.createGoMatrix()
  
  LFC_table <- data.frame(WT1_4w_LFC=WT1_4w$log2FoldChange,WT1_12w_LFC=WT1_12w$log2FoldChange,row.names = rownames(WT1_4w))
  x<- merge(LFC_table,geneID2entrez,by.x = 0, by.y = "ensembl_gene_id",all = F)
  x <- x[!is.na(x$entrezgene),]
  x <- ddply(x[,-1],"entrezgene",numcolwise(sum))
  LFC_tableENT <- data.frame(early_FSGS=x$WT1_4w_LFC,late_FSGS=x$WT1_12w_LFC,row.names = x$entrezgene)
  
  lfc <- t(sapply(RobertGO_sigUni, function(x) {set <- gomatrix[which(rownames(gomatrix)==x),]
  set <- names(set[set==1])
  colMeans(LFC_tableENT[rownames(LFC_tableENT)%in%set,])}))
  
  lfc <- cbind.data.frame(lfc,union(RobertGO_4w$Term[which(RobertGO_4w$Primary.Fisher.adj<0.1)],RobertGO_12w$Term[which(RobertGO_12w$Primary.Fisher.adj<0.1)]))
  colnames(lfc) <- c("meanLFC_earlyFSGS","meanLFC_lateFSGS","Primary")
  write.table(lfc,sep = "\t",quote = F,file="/data/user/tpadvits/PROJECTS/MARTIN_PJ/RNAseq_analysis/WT1/DEexon_intron/LFCheatmap_GO/robertGO_clustAll_LFCcompare.csv")
}

####======== Signaling Pathway Impact Analysi (SPIA) ======  
#======================================================#
############ SPIA via graphite
### prepare data and run SPIA
{ 
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "mmusculus_gene_ensembl")
  geneID2entrez=biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id"), mart = mart)
  library("graphite")
  library("SPIA")
  ## download prepared KEGG and REACTOME pathways to run with graphite
  path_REACT <-pathways("mmusculus", "reactome")
  path_REACT <- convertIdentifiers(path_REACT, "ENTREZID")
  path_KEGG <- pathways("mmusculus", "kegg")
  # prepare for SPIA
  prepareSPIA( path_KEGG, "KEGGAll")
  prepareSPIA( path_REACT, "reactomeAll")
  
  
  ## databases from 10.2017
  ## databases from 27-04-2018
  
  ## define gene set
  # results for 4w, 12w and aging have shrunk LFC, disease progression is with non-shrunk lfc due to the use of contrasts
  # updates 02.05.19
  WT1_4w=read.table( "/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/DE_analysis/DEexonIntron_4w_glvl.csv" )
  WT1_12w=read.table( "/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/DE_analysis/DEexonIntron_12w_glvl.csv" )
  # WT1_disprog=read.table(row.names = 1 , header = T, "/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/DE_analysis/DEexonIntron_dprogOnly_glvl.csv")
  # WT1_aging=read.table(row.names = 1 , header = T,"/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/DE_analysis/DEexonIntron_aging_glvl.csv")
  WT1_4w$gName <- WT1_12w$gName <- WT1_aging$gName <- WT1_disprog$gName <- tx2gene$external_gene_name[ match( rownames(WT1_aging) , tx2gene$ensembl_gene_id ) ]
  Nphs2_4w <- readRDS( file="/data/user/tpadvits/PROJECTS/MARTIN_PJ/RNAseq_analysis/NPHS2/Nphs2_IntExon_glvl_4wDE.rda" )
  Nphs2_12w <- readRDS( file="/data/user/tpadvits/PROJECTS/MARTIN_PJ/RNAseq_analysis/NPHS2/Nphs2_IntExon_glvl_12wDE.rda" )

  inputSPIA <- list( Nphs2_4w, Nphs2_12w )
  inputSPIA <- lapply( inputSPIA, as.data.frame )
  names(inputSPIA) <- c( "Nphs2_4w" , "Nphs2_12w" )
    

  SPIAonList <- function(inputSPIA , strict = F ) 
    { 
    ensGene_entres <- read.table(file='/data/user/tpadvits/PROJECTS/MARTIN_PJ/RNAseq_analysis/WT1/kallisto_reverse/cDNAgenes_entrezID.csv',sep="\t",header = T)
    
    x <- inputSPIA
      
      if(strict!=T){
        ## define background gene set
        ALL_genes <- unique( ensGene_entres$Entrez.Gene.ID )
        ALL_genes <- union( ALL_genes, unique( geneID2entrez$entrezgene[ which( geneID2entrez$ensembl_gene_id%in%rownames(x) |  geneID2entrez$external_gene_name%in%rownames(x) ) ]) )
        ALL_genes <- paste("ENTREZID", ALL_genes, sep = ":") # modify gene names 
        
      } else if (strict==T){
        ## OR define background gene set more strictly as only genes participated in DE analysis
        ALL_genes <- unique( geneID2entrez$entrezgene[ which(geneID2entrez$ensembl_gene_id%in%rownames(x)) ] )
        ALL_genes <- paste("ENTREZID", ALL_genes, sep = ":") # modify gene names 
      }
      
      ## DE at week 4, 12, disProg or aging in entrez IDs
      # modify gene names for compatibility with GRAPHITE
      x <- x[which(x$padj<0.05),]
      x$entrezgene <- geneID2entrez$entrezgene_id[match( rownames(x), geneID2entrez$ensembl_gene_id)]
      x <- aggregate( log2FoldChange~entrezgene,x,sum)
      x <- setNames( x$log2FoldChange,paste("ENTREZID", x$entrezgene, sep = ":"))
      # # make sure no genes are missing from the background
      # ALL_genes <- union(ALL_genes, names(x) )
      
      ### run SPIA
      SPIA_REACT <- runSPIA(de=x, all=ALL_genes, pathwaySetName="reactomeAll",verbose=F,nB=2000)
      # SPIA_REACT_sig=SPIA_REACT[which(SPIA_REACT$pGFWER<=0.05),]
      SPIA_KEGG <- runSPIA(de=x, all=ALL_genes, pathwaySetName="KEGGAll",verbose=F,nB=2000)
      # SPIA_KEGG_sig=SPIA_KEGG[which(SPIA_KEGG$pGFWER<=0.05),]
      
      return(list(SPIA_REACT,SPIA_KEGG))
    }


  SPIA_both_NONstrictBkg <- SPIAonList(strict = F ) 
  #SPIA_both_strictBkg <- SPIAonList(strict = T ) 
  
  # saveRDS(SPIA_both_NONstrictBkg, file="/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/SPIA/allLFCshr_SPIA_both.NONstrictBkg.rda")
  # saveRDS(SPIA_both_strictBkg, file="/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/SPIA/allLFCshr_SPIA_both.strictBkg.rda")
  
  # SPIA_REACTnew <- readRDS( file="/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/SPIA/allLFCshr_SPIA_both.NONstrictBkg.rda")[[1]]
  # SPIA_KEGGnew <- readRDS( file="/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/SPIA/allLFCshr_SPIA_both.NONstrictBkg.rda")[[2]]
  
  # all_SPIA_REACTallbkg <- readRDS(  "/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/SPIA/all_SPIA_REACT.allannotbkg.rda" )
  # all_SPIA_KEGGallbkg <- readRDS(  "/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/SPIA/all_SPIA_KEGG.allannotbkg.rda" )
  
}

### Venn diagram to compare 4 VS 12 weeks pathways
{ library(VennDiagram)
  test_compare_gene0.1 <- list(early_FSGS=SPIA_KEGG_w4$Name[which(SPIA_KEGG_w4$pGFWER<0.1)],
                               late_FSGS=SPIA_KEGG_w12$Name[which(SPIA_KEGG_w12$pGFWER<0.1)])
  venn.plot <- venn.diagram(test_compare_gene0.1 , NULL, fill=c("darkblue","gray"), alpha=c(0.5,0.5), 
                            cex = 4, cat.fontface=4, category.names=c("early_FSGS\n N=34", "late_FSGS\n N=26"), main="overlap of KEGG pathways\n identified by SPIA (pGFWER<0.01)")
  # To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
  grid.draw(venn.plot)
}

### wrapper function to plot Heatmap results of SPIA 
path_sel <- read.table(sep = "\t","/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/SPIA/SPIAKEGGREACTOMEOVERLAP_2018.csv")
path_sel2019 <- path_sel[-c(26,31),]
# combine REACTOME and KEGG results and add dbID to pathway names
set4w_both <- rbind(SPIA_REACTlist[[1]],SPIA_KEGGnew[[1]])
set4w_both$Name <- paste(set4w_both$Name,c(rep("REACTOME", dim(SPIA_REACTlist[[1]])[1]),rep("KEGG",dim(SPIA_KEGGnew[[1]])[1])),sep = "^")
set12w_both <- rbind(SPIA_REACTlist[[2]],SPIA_KEGGnew[[2]])
set12w_both$Name <- paste(set12w_both$Name,c(rep("REACTOME", dim(SPIA_REACTlist[[2]])[1]),rep("KEGG",dim(SPIA_KEGGnew[[2]])[1])),sep = "^")

# load SPIA results 
ll <- list.files(pattern = "^SPIA_" , "/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/SPIA",full.names = T)
SPIA_results <- lapply(ll, readRDS)
SPIA_REACT <-  lapply(SPIA_results, `[[`, 1)
SPIA_REACT[c(1,2)] <- SPIA_REACT[c(2,1)]
SPIA_KEGG <-  lapply(SPIA_results, `[[`, 2)
SPIA_KEGG[c(1,2)] <- SPIA_KEGG[c(2,1)]

# SPIA_REACT <- list(all_SPIA_REACTallbkg[[1]],all_SPIA_REACTallbkg[[2]],SPIA_REACTnew[[3]],SPIA_REACTnew[[4]])
# SPIA_KEGG <- list(all_SPIA_KEGGallbkg[[1]],all_SPIA_KEGGallbkg[[2]],SPIA_KEGGnew[[3]],SPIA_KEGGnew[[4]])

# separate heatmaps 
# path_sel1 <- path_sel[which(paste(path_sel[,1],path_sel[,2],sep = "^")%in%union_tab$union[which(union_tab[,2]<=0.1 & union_tab[,4]<=0.1)]),]
# path_sel2 <- path_sel[which(paste(path_sel[,1],path_sel[,2],sep = "^")%in%union_tab$union[which(union_tab[,2]<0.1 & union_tab[,4]>0.1)]),]
# path_sel3 <- path_sel[which(paste(path_sel[,1],path_sel[,2],sep = "^")%in%union_tab$union[which(union_tab[,2]>0.1 & union_tab[,4]<0.1)]),]

nameCompare <-  c("FSGS_4w","FSGS_12w","aging","dis.prog.","dis.prog.(-aging)")

plotSPIA <- function( set1, set2, SPIA_KEGGlist=NULL, SPIA_REACTlist=NULL, nameL = 50 , DendroBar=F, ringPlot=F,  scaleR=c(-10,10),
                      alpha=0.1, path_type="path_type", selected=NULL, C_height = 5 , orderLFC=F , colab=nameCompare)
 { 
  require(ggplot2)
  require(reshape)
  require(plyr)
  
  if (is.null(SPIA_KEGGlist)){

   if (is.null(selected)){
    set1_sig <- set1[which(set1$pGFWER<alpha),]
    set2_sig <- set2[which(set2$pGFWER<alpha),]
  } else {
    selected <- paste(selected[,1],selected[,2],sep = "^")
    names(set1) <- names(set2) 
    set1_sig <- set1[which(set1[,1]%in%selected),] 
    set2_sig <- set2[which(set2[,1]%in%selected),] 

  }
  
  union <- union(set1_sig$Name, set2_sig$Name)
  
  union_tab <- cbind(union,set1[(match(union,set1$Name)),c(9,10)],set2[(match(union,set2$Name)),c(9,10)])
  union_tab[is.na(union_tab[,2]),2] <- 1  # eleminate possible NAs
  union_tab[is.na(union_tab[,4]),4] <- 1
  union_tab[is.na(union_tab[,3]),3] <- "Inhibited"
  union_tab[is.na(union_tab[,5]),5] <- "Inhibited"
  names(union_tab)=c(path_type,"early_FSGS_pGFWER","set1_status","late_FSGS_pGFWER","set2_status")
  
  X <- union_tab
  names <- gsub("_pGFWER","",names(X)[c(1,2,4)])
  
  x <- revalue(as.factor(X[,3]), c("Activated"=-1, "Inhibited"=1))
  x =as.numeric(levels(x))[x]
  y <- revalue(as.factor(X[,5]), c("Activated"=-1, "Inhibited"=1))
  y= as.numeric(levels(y))[y]
  
  inter_toPlot=data.frame(both=X[,1],FSGS_4w=log10(X[,2])*x, FSGS_12w=log10(X[,4])*y)
  } else { 
    if(is.null(selected)) {stop("provide a list of selected patways to plot SPIA results for all comparisons")
      } else{
      # ll <- length(SPIA_KEGGlist[[1]][match( selected[,1] , SPIA_KEGGlist[[1]]$Name),1] )
      path_union <-   paste(selected[,1],selected[,2],sep = "^")
      
      # check if length of Reactome and SPIA are equal if no throw an error
      if(length(SPIA_KEGGlist)==length(SPIA_REACTlist)) {
        
      setAll_sig <- lapply(seq(SPIA_KEGGlist) , function(x){
          set_sig <- rbind( SPIA_KEGGlist[[x]] , SPIA_REACTlist[[x]] )
          set_sig$Name <- c( paste( set_sig$Name[1:dim(SPIA_KEGGlist[[x]])[1]], '^KEGG' , sep = "") , paste(set_sig$Name[(dim(SPIA_KEGGlist[[x]])[1]+1):dim(set_sig)[1]], '^REACTOME' , sep = "") )
          set_sig <- set_sig[ match( path_union , set_sig$Name),] 
          return(set_sig)
        })
        

        ### union of all comparisons
        union_tab <- lapply(setAll_sig , function(x) x[(match(path_union,x[,1])),c(9,10)])
        union_tab <- Reduce( cbind , union_tab) 
        union_tab <- cbind.data.frame( path_union, union_tab )
        
        union_tab[ , seq( 2 ,(dim(union_tab)[2]-1), by=2)][is.na(union_tab[ , seq( 2 ,(dim(union_tab)[2]-1), by=2)] )] <- 1  # eleminate possible NAs
        union_tab[ , seq( 3 ,(dim(union_tab)[2]), by=2)][is.na(union_tab[ , seq( 3 ,(dim(union_tab)[2]), by=2) ] )] <- "Inhibited"
        
        X <- union_tab
        # revalue in a list
        X[ , seq( 3 , (dim(X)[2]) , by=2) ] <- lapply( X[, seq( 3 , (dim(X)[2]) , by=2)] , function(y) {
          x <- revalue(as.factor(y), c("Activated"=-1, "Inhibited"=1)) 
          x <- as.numeric(levels(x))[x]
          return(x) })
        
        inter_toPlot <- cbind.data.frame( X[,1] , log10(X[ , seq( 2 ,(dim(union_tab)[2]-1), by=2)] ) * X[, seq( 3 ,(dim(union_tab)[2]), by=2)] )
        names(inter_toPlot) <- c( "both" , colab )
      } else stop("lists of SPIA results for KEGG and REACTOME should be of the same length")
     
    }
    
  }

    # round significance values (in case we want to display them)
    inter_toPlot[,-1] <- round((inter_toPlot[,-1]),2)
    inter_toPlot[,1] <- sapply(inter_toPlot[,1], function(x) wrap_text(as.character(x),n=nameL) )
    
    if (!is.null(scaleR)) {
      if (length(scaleR)==2 & class(scaleR)=="numeric"){
        #inter_toPlot <- cbind.data.frame(both=inter_toPlot[,1],t(scale(t(inter_toPlot[,-1])) ))
        inter_toPlot[,-1][inter_toPlot[,-1] < (scaleR[1])] <- scaleR[1]
        inter_toPlot[,-1][inter_toPlot[,-1] > (scaleR[2])] <- scaleR[2]
      } else stop("scaleR argument is incorrect: set to NULL or supply a numeric vector of length 2 with lower and upper limits for lfc scale")
        
  }
   
    if(orderLFC==F){
      ordR <- hclust( dist(inter_toPlot[,-1] , method = "euclidean"), method = "ward.D" )$order
    } else  ordR <- order(rowMeans(inter_toPlot[,-1]))
    

    # cluster pathways
    # melt for use with ggplots
    melted_cormat <- melt(inter_toPlot, na.rm = TRUE)
    melted_cormat[,1] <- factor(melted_cormat[,1], levels = unique(melted_cormat[,1]))
    melted_cormat[,4] <- abs( melted_cormat[,3])

    # order rows
    melted_cormat[,1] <- factor( melted_cormat[,1] , levels = inter_toPlot[,1][ordR], labels =   inter_toPlot[,1][ordR] )

    if(ringPlot==F) {
      
      # if(value==T) # use of values in tiles is depricated for now due to no utility
      #   p1 <- ggplot(data = melted_cormat, aes(variable, melted_cormat[,1], fill = value)) + geom_tile(color = "white",stat = "identity",size=0)+
      #   geom_text(aes(label=V4),cex = 2.5) +  scale_y_discrete(expand=c(0,0) , position = "right") + # to show values in tiles
      #   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
      #                        midpoint = 0, limit = c(min(melted_cormat$value),max(melted_cormat$value),labels=NULL), space = "Lab", 
      #                        name="Direction of pathway change\n pseudo-colours") +
      #   theme_minimal()+  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle = 45,hjust=1, size=14),
      #                           axis.title.y=element_blank())  + coord_fixed(ratio=0.2)
      # else
      
      p1 <- ggplot(data = melted_cormat, aes(variable, melted_cormat[,1], fill = value)) + geom_tile(color = "white",stat = "identity",size=0)+
        scale_y_discrete(expand=c(0,0) , position = "right") + scale_x_discrete(expand=c(0,0) ) +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                             midpoint = 0, labels=abs , space = "Lab", 
                             limit = c(min(melted_cormat$value),max(melted_cormat$value)),
                             name="Pathway change\nlog10(pGFWER)\nin pseudo-colours") +
        theme_minimal()+  theme(axis.title=element_blank(), 
                                element_line(colour=NULL) ,
                                axis.ticks= element_blank(),
                                axis.text.x=element_text(angle = 90, hjust=1 , vjust=0.5 , size=14 ),
                                axis.text.y = element_text(size=10 ), axis.title.y=element_blank())  + coord_fixed(ratio=0.2)
      
      if(DendroBar==T) {
        
        # make cluster bar in addition to a heatmap
        require(ggdendro)
        require(gridExtra)
        require(ggpubr)
        
      hc <- hclust( dist(as.data.frame( inter_toPlot[,-1] , row.names = as.character( inter_toPlot[,1]) ) , method = "euclidean"), method = "ward.D" )
      df2 <- data.frame( cluster=cutree(hc,C_height) , states=factor(hc$labels,levels=hc$labels[hc$order]))
      # plot
      p2 <- ggplot( df2, aes(x=1,y=states,fill=factor(cluster))) + geom_tile()+
        scale_y_discrete(expand=c(0,0))+
        scale_x_discrete(expand=c(0,0),name = NULL ) +
        theme( axis.title.y = element_blank() , element_line(colour="white") ,
               axis.ticks= element_blank(),
               axis.text= element_blank(),
               axis.title.x = element_text(face="bold" , angle = 90 , size=14 ,  vjust = 0.5),
               legend.position="none" )
      ## arrange grids
      ggpubr::ggarrange( p2, p1, ncol=2 , widths=c(1,20) , align = "h")
      
     } else(p1)
    } else {
      require(forcats)
      sequence_length = length(unique(melted_cormat[,1]))
      first_sequence = c(1:(sequence_length%/%2)) 
      second_sequence = c((sequence_length%/%2+1):sequence_length) 
      first_angles = c(90 - 180/length(first_sequence) * first_sequence)
      second_angles = c(-90 - 180/length(second_sequence) * second_sequence)
      
      melted_cormat$variable <- fct_rev(melted_cormat$variable)
        
      ggplot2::ggplot(data = melted_cormat, aes(x=variable, y=melted_cormat[,1], fill = value )) + geom_tile(color = "white",stat = "identity",size=0)+
        scale_x_discrete(expand=c(1,1)) +  
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                             midpoint = 0, labels=abs , space = "Lab", 
                             #limit = c(min(melted_cormat$value),max(melted_cormat$value)),
                             name="Pathway change\nlog10(pGFWER)\nin pseudo-colours") + 
         theme_minimal()+  
          theme(axis.title=element_blank(), plot.margin=unit(c(5,5,5,5),"cm") ,
                panel.border = element_blank(),
                #panel.grid.major = element_blank(),
                                element_line(colour=NULL) ,
                                axis.ticks= element_blank(),
                                axis.text.x = element_text(angle = c(first_angles,second_angles), hjust=rep(1,length(melted_cormat[,1])) , vjust=-1 , size=10 ),
                                axis.text.y = element_blank(), axis.title.y=element_blank()) +  
        coord_polar(theta="y", clip="off") # circularize plot, dont clip pathway names: clip="off"
        
        }
}


# makle a plot (circular)
plotSPIA(   set1= set4w_both, set2=set12w_both,  DendroBar=F, ringPlot=F, scaleR=NULL, path_type="path_type", 
            alpha=0.1 , C_height = 5, nameL = 200 ,  orderLFC=F )
# non circular
plotSPIA(   SPIA_KEGGlist=SPIA_KEGG , SPIA_REACTlist=SPIA_REACT ,  DendroBar=F, ringPlot=F, scaleR=c(-15,15) , path_type="path_type", 
            selected = path_sel2019, C_height = 5, nameL = 200 ,  orderLFC=F )

### heatmap for selected paths
 {    
# SPIA_4w <- rbind(SPIA_KEGG_w4,SPIA_REACT_w4)
#   SPIA_4w[,11] <- c(rep("KEGG",dim(SPIA_KEGG_w4)[1]),rep("REACT",dim(SPIA_REACT_w4)[1]))
#   SPIA_4w <- subset(SPIA_4w, (SPIA_4w$Name%in%path_sel$V1[which(path_sel$V2=="KEGG")] & SPIA_4w$V11=="KEGG") | (SPIA_4w$Name%in%path_sel$V1[which(path_sel$V2=="Reactome")] & SPIA_4w$V11=="REACT"))
#   SPIA_4w$Name[SPIA_4w$V11=="REACT"] <- paste(SPIA_4w$Name[SPIA_4w$V11=="REACT"],"*",sep = "")
#   SPIA_12w <- rbind(SPIA_KEGG_w12,SPIA_REACT_w12)
#   SPIA_12w[,11] <- c(rep("KEGG",dim(SPIA_KEGG_w12)[1]),rep("REACT",dim(SPIA_REACT_w12)[1]))
#   SPIA_12w <- subset(SPIA_12w, (SPIA_12w$Name%in%path_sel$V1[which(path_sel$V2=="KEGG")] & SPIA_12w$V11=="KEGG") | (SPIA_12w$Name%in%path_sel$V1[which(path_sel$V2=="Reactome")] & SPIA_12w$V11=="REACT"))
#   SPIA_12w$Name[SPIA_12w$V11=="REACT"] <- paste(SPIA_12w$Name[SPIA_12w$V11=="REACT"],"*",sep = "")
# 
#   # prepare data for plotting
#   library(plyr)
#    union <- union(SPIA_4w$Name,SPIA_12w$Name)
#   union_tab <- cbind(union,SPIA_4w[(match(union,SPIA_4w$Name)),c(9,10)],SPIA_12w[(match(union,SPIA_12w$Name)),c(9,10)])
#   union_tab[is.na(union_tab[,2]),2] <- 1  # eleminate possible NAs
#   union_tab[is.na(union_tab[,4]),4] <- 1
#   union_tab[is.na(union_tab[,c(3)]),c(3)] <- "Inhibited"
#   union_tab[is.na(union_tab[,c(5)]),c(5)] <- "Inhibited"
#   names(union_tab)=c("pathway","early_FSGS_pGFWER","set1_status","late_FSGS_pGFWER","set2_status")
#   union_tab[union_tab$early_FSGS_pGFWER<0.1 & union_tab$late_FSGS_pGFWER<0.1,6] <- "both"
#   union_tab[union_tab$early_FSGS_pGFWER<0.1 & union_tab$late_FSGS_pGFWER>=0.1,6] <- "early_only"
#   union_tab[union_tab$early_FSGS_pGFWER>=0.1 & union_tab$late_FSGS_pGFWER<0.1,6] <- "late_only"
#   
#   #  x <- union_tab[union_tab$V6=="both",]
#   # names <- gsub("_pGFWER","",names(x)[c(1,2,4)])
#   # inter_toPlot=x[order(x[,3],x[,2]),]     # order according to status and p-value of the first condition
#   # x=revalue(as.factor(inter_toPlot[,3]), c("Activated"=-1, "Inhibited"=1))
#   # x=as.numeric(levels(x))[x]
#   # y=revalue(as.factor(inter_toPlot[,5]), c("Activated"=-1, "Inhibited"=1))
#   # y=as.numeric(levels(y))[y]
#   # inter_toPlot=data.frame(inter_toPlot[,1],log10(inter_toPlot[,2])*x,log10(inter_toPlot[,4])*y)
#   # names(inter_toPlot) <- names
#   # inter_toPlot=inter_toPlot[order(inter_toPlot[,2],decreasing = T),]
#   # inter_toPlot[,2:3] <- round((inter_toPlot[,2:3]),2)
#   
#   union_tab_trans <- 
#   
#   melted_cormat <- melt( union_tab[,c(1,2,4)] , na.rm = TRUE)
#   # melted_cormat[,1] <- factor(melted_cormat[,1], levels = unique(melted_cormat[,1]))
#   # melted_cormat[,4] <- abs(melted_cormat[,3])
#   melted_cormat$variable <- factor(melted_cormat$variable, labels = c("early FSGS","late FSGS"))
#   
#   
#   labl <- c(20,"10 down-regulated",5,0,5,"10 up-regulated")
#   
#   ordR <- hclust( dist(XX , method = "euclidean"), method = "ward.D" )$order
#   ordC <- hclust( dist( t(XX) , method = "euclidean"), method = "ward.D" )$order
#   
#   XX.m <- melt( as.matrix( XX ) )
#   # # add significance marks to tiles
#   XX.m$padj <- melt(HomoMouse_qval)$value
#   XX.m$padj <-  ifelse( is.na(XX.m$padj) , "", ifelse( XX.m$padj > 0.1 , "" , ifelse( XX.m$padj > 0.01 ,  "*" , "**") ) )
#   
#   # order rows
#   XX.m$X1 <- factor( XX.m$X1, levels = rownames(XX)[ordR], labels =  rownames(XX)[ordR] )
#   # order columns
#   XX.m$X2 <- factor( XX.m$X2, levels = colnames(XX)[ordC], labels =  colnames(XX)[ordC] )
#   
#   g1<- ggplot(data = melted_cormat, aes(variable, melted_cormat[,1], fill = value)) + geom_tile(color = "white",stat = "identity",size=0)+
#     scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                          midpoint = 0, limit = c(min(melted_cormat$value),max(melted_cormat$value)),labels=labl, space = "Lab", 
#                          name="-log10*q-value\n of pathway change\n in pseudo-colours") +
#     theme_minimal()+  theme(text = element_text(size=14),axis.title.x=element_blank(),axis.text.x=element_text(angle = 45,hjust=1),
#                             axis.title.y=element_blank())  + coord_fixed(ratio=0.2)
#   # p1 <- as.grob(g1)
#   # p2 <-  as.grob(g2)
#   # p3 <-  as.grob(g3)
#   # grid.arrange(grobs=list(p1, p2, p3),ncol=1,heights = c(2, 1, 1))
 } 
   
### venn of pathway overlap
  library(VennDiagram)
    test_compare_gene0.1 <- list(early_FSGS=SPIA_4w$Name[which(SPIA_4w$pGFWER<0.1)],
                                 late_FSGS=SPIA_12w$Name[which(SPIA_12w$pGFWER<0.1)])
    venn.plot <- venn.diagram(test_compare_gene0.1 , NULL, fill=c("green","salmon"), alpha=c(0.5,0.5),
                              cex =4,cat.fontface=4, cat.cex=2,cat.just=list(c(0.2,-0.35),c(0.8,-0.35)),category.names=c("early FSGS\n N=26", "late FSGS\n N=25"))
    # To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
    grid.draw(venn.plot)


### pathview visualization of KEGG pathways highlited by SPIA ###
  { 
    
    library(plyr)
  # make an LFC table with ENTREZ ids
  LFC_table <- data.frame(WT1_4w_LFC=WT1_4w$log2FoldChange,WT1_12w_LFC=WT1_12w$log2FoldChange,row.names = rownames(WT1_4w))
  x<- merge(LFC_table,geneID2entrez,by.x = 0, by.y = "ensembl_gene_id",all = F)
  x <- x[!is.na(x$entrezgene),]
  x <- ddply(x[,-1],"entrezgene",numcolwise(sum))
  LFC_tableKEGG <- data.frame(WT1_4w_LFC=x$WT1_4w_LFC,WT1_12w_LFC=x$WT1_12w_LFC,row.names = x$entrezgene)
  # select KEGG pathways to vizualize
  set1_sig <- set1[which(set1$pGFWER<0.1),]
  set2_sig <- set2[which(set2$pGFWER<0.1),]
  union <- union(set1_sig$Name,set2_sig$Name)
  
  setwd("/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/SPIA/KEGG_PATHVIEW/")  
  for(i in 1:length(union)){
    kegg_path <- path_KEGG[[which(names(path_KEGG)==union[i])]]
    
    # pathview(gene.data=LFC_tableKEGG, pathway.id=sub("mmu:","",kegg_path@id), species="mmu",kegg.native=T, out.suffix = paste(kegg_path@title,"FSGS",sep="."),
    #          kegg.dir = "/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/RNAseq/SPIA/PATHVIEW/KEGG_path",bins=list(gene=20, cpd=10))
    pathview(gene.data=LFC_tableKEGG, pathway.id=sub("mmu:","",kegg_path@id), species="mmu",kegg.native=F, out.suffix = paste(kegg_path@title,"FSGS",sep="."),
             kegg.dir = "/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/RNAseq/SPIA/KEGG_PATHVIEW/KEGG_path",bins=list(gene=20, cpd=10))
  }
}  

### 2-way evidence SPIA plot
      XX <- "/data/user/tpadvits/PROJECTS/MARTIN_PJ/sc_RNAseq/snRNAseq_kif3a/kallisto/SPIA/" 
      spia2way_ggplot <- function( pathSPIA = XX , namesSPIA=nameCompare , pathtype="REACTOME",  SPIA_list = SPIA_REACT, combinemethod = "fisher", threshold = 0.1) 
        {
        library(SPIA)
        library(ggplot2)
        library(ggrepel)
        # helper function from SPIA
        getP2<-function(pG, combine=combinemethod)
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
        
        for (i in 1:length(SPIA_list))
        {
          
          x <- SPIA_list[[i]]
          x$pPERT[is.na(x$pPERT)] <- 1
          x$pNDE[is.na(x$pNDE)] <- 1
          
          
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
            scale_colour_manual(name = 'Pathway status', values = setNames(c("dark grey","blue","red") ,c("dark grey","blue","red")), labels=c('down-regulated','non-sig','up-regulated'))+
            geom_abline(aes(intercept = -log(getP2(tr1,combinemethod)^2), slope = -1), color = "green" ,lwd=1)+
            geom_abline(aes(intercept = -log(getP2(tr2,combinemethod)^2), slope = -1), color = "grey50" ,lwd=1) + 
            # ylim(0, 15) + xlim(0, 15) + 
            theme_minimal() + theme(plot.title = element_text(face="bold",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"), 
                                    legend.text = element_text(size = 20),legend.title = element_text(size=20),legend.key.size = unit(2, "cm"),strip.text.y = element_text(size = 20))+
            guides(colour = guide_legend(override.aes = list(size=20))) + # control size of of points in legend 
            geom_text_repel( size = 6 , point.padding =1 , aes(label = x$NameSig),box.padding   = 1, segment.color = 'grey50')  #ensure that labels are not jammed together
          
          pdf(file = paste(pathSPIA,namesSPIA[i],"_SPIAevid_",pathtype,".pdf",sep = ""), width = 15, height = 12) 
          print (PP)
          dev.off()
          
          png(file = paste(pathSPIA,namesSPIA[i],"_SPIAevid_",pathtype,".png",sep = ""), width = 1000, height = 700) 
          print (PP)
          dev.off()
        }
      }
     

### compare 2 conditions using SPIA
    {
      
      SPIA_Nphs2 <- readRDS("/data/user/tpadvits/PROJECTS/MARTIN_PJ/RNAseq_analysis/NPHS2/SPIA/Nphs2_SPIA.nonstrict.rda")
      # path_sel <-  read.table(sep = "\t","SPIAKEGGREACTOMEOVERLAP_2018.csv")
      # path_sel <- paste(path_sel[,1],path_sel[,2],sep = "^")
      
    # define string wrapper function
      wrap_text <- function(string, n=40) 
      { require(stringr)
        if(nchar(string) > n){
          spaces <- str_locate_all(string, " ")[[1]][,1]
          chars  <- nchar(string)
          for(i in 1:floor(chars/n)) {
            s <- spaces[which.min(abs(spaces - n*i))]
            substring(string, s, s) <- "\n "
          }
          return(string)
        } else string <- string
        
      }
      
      plotSPIA_2D <- function( SPIA_results=SPIA_Nphs2 , nameL = 50 ,
                            alpha=0.1 , selected=NULL , Nsiglabel=20 , plotTitle="")
      {    
        #path_sel <- path_sel[-c(26,31),]
        set4w_both <- rbind( SPIA_results[[2]][[1]] , SPIA_results[[1]][[1]] ) 
        set4w_both$Name <- paste( set4w_both$Name, c(rep("REACTOME", dim( SPIA_results[[2]][[1]])[1] ),rep("KEGG",dim(SPIA_results[[1]][[1]])[1] )),sep = "^") 
        set12w_both <- rbind( SPIA_results[[2]][[2]] , SPIA_results[[1]][[2]] ) 
        set12w_both$Name <- paste(set12w_both$Name,c(rep("REACTOME", dim( SPIA_results[[2]][[2]])[1] ),rep("KEGG",dim(SPIA_results[[1]][[2]])[1] )),sep = "^")
        
      set1_sig <- set4w_both[which(set4w_both$pGFWER<alpha),]
      set2_sig <- set12w_both[which(set12w_both$pGFWER<alpha),]
      
      union <- union(set1_sig$Name, set2_sig$Name)
      
      union_tab <- cbind( union, set4w_both[(match(union,set4w_both$Name)),c(9,10)],set12w_both[(match(union,set12w_both$Name)),c(9,10)])
      union_tab[is.na(union_tab[,2]),2] <- 1  # eleminate possible NAs
      union_tab[is.na(union_tab[,4]),4] <- 1
      union_tab[is.na(union_tab[,3]),3] <- "Inhibited"
      union_tab[is.na(union_tab[,5]),5] <- "Inhibited"
      
      X <- union_tab
      x <- revalue(as.factor(X[,3]), c("Activated"=-1, "Inhibited"=1))
      x <- as.numeric(levels(x))[x]
      y <- revalue(as.factor(X[,5]), c("Activated"=-1, "Inhibited"=1))
      y <- as.numeric(levels(y))[y]
      
      inter_toPlot=data.frame( both=X[,1],FSGS_4w=log10(X[,2])*x, FSGS_12w=log10(X[,4])*y)
      inter_toPlot$shapeD <- as.factor( ifelse( abs(inter_toPlot$FSGS_4w) > 1 & abs(inter_toPlot$FSGS_12w) > 1, 'change in both' ,ifelse(  abs(inter_toPlot$FSGS_4w) < 1 & abs(inter_toPlot$FSGS_12w) > 1 , '12w only' , '4w only' ) ) )
      

      inter_toPlot$colorD <- factor(  inter_toPlot$shapeD , levels = c('change in both','4w only', '12w only') , labels = c('green','magenta','orange'))
   
      
      inter_toPlot$both <- ifelse( inter_toPlot$both%in% union( set1_sig$Name[order(set1_sig$pGFWER)][ seq(Nsiglabel) ] , set2_sig$Name[order(set2_sig$pGFWER)][ seq(Nsiglabel) ]) , as.character(inter_toPlot$both) , "")
      
      # wrap long pathway names
      inter_toPlot[,1] <- sapply(inter_toPlot[,1], function(x) wrap_text(as.character(x),n=nameL) )
      
      
      ggplot(data = inter_toPlot, aes(x=FSGS_4w, y=FSGS_12w, shape=shapeD)) + geom_hline( yintercept=0 ) + geom_vline( xintercept=0 )  +
        geom_point(size=5,  fill=inter_toPlot$colorD  ) + 
        scale_shape_manual("", values = c( 'change in both'=21,'4w only'=22, '12w only'=24)) + ggtitle(plotTitle) +
        geom_text_repel(data=inter_toPlot,aes(label=both), box.padding   = 0.35, point.padding = 0, segment.color = 'grey50', size=4)  + #ensure that labels are not jammed together
        theme_minimal() +  theme(axis.text=element_text(size=20), plot.title = element_text(size = 24, face = "bold", hjust = 0.5) ,
                                 axis.title=element_text( size=24,face="bold"), legend.text = element_text(size = 24),legend.title = element_text(size=24) )+
        scale_size(range = c(5, 6)) + theme(legend.key.size = unit(2.5, "cm")) + labs(colour = "condition") # +
      }
      
      
    }
    
####======== ENRICHr - web-based tool ===================
## http://amp.pharm.mssm.edu/Enrichr/



####======== 2D plots ==========
######### 2D plot for GO results
robertGO_clustAll_compare <- read.table(header = T,sep = "\t", "/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/LFCheatmap_GO/robertGO_clustAll_LFCcompare.csv",stringsAsFactors = F )

temp_var <- predict.lm(lm(robertGO_clustAll_compare$meanLFC_lateFSGS~robertGO_clustAll_compare$meanLFC_earlyFSGS), 
                       interval="prediction" , level = 0.95 )
resi <- resid(lm(robertGO_clustAll_compare$meanLFC_lateFSGS~robertGO_clustAll_compare$meanLFC_earlyFSGS))  

new_df <- cbind.data.frame(robertGO_clustAll_compare, temp_var,resi)
new_df$color=1
new_df$color[new_df$meanLFC_lateFSGS<new_df$lwr & new_df$meanLFC_lateFSGS<0 & new_df$meanLFC_earlyFSGS<0]<-2
new_df$color[new_df$meanLFC_lateFSGS>new_df$upr & new_df$meanLFC_lateFSGS<0 & new_df$meanLFC_earlyFSGS<0]<-3
new_df$color[new_df$meanLFC_lateFSGS>new_df$upr & new_df$meanLFC_lateFSGS>0 & new_df$meanLFC_earlyFSGS>0]<-4
new_df$color[new_df$meanLFC_lateFSGS<new_df$lwr & new_df$meanLFC_lateFSGS>0 & new_df$meanLFC_earlyFSGS>0]<-5
new_df$color[((new_df$meanLFC_lateFSGS<new_df$lwr & new_df$meanLFC_lateFSGS<0 & new_df$meanLFC_earlyFSGS>0) | (new_df$meanLFC_lateFSGS>new_df$upr & new_df$meanLFC_lateFSGS>0 & new_df$meanLFC_earlyFSGS<0))]<-6

revalue(x, c("beta"="two", "gamma"="three"))

library(plyr)

new_df$status <- revalue( as.factor(new_df$color) , c("1"= "similarily changing at both stages","2"= "more downreg. in late FSGS" ,"3"= "more downreg. in early FSGS","4"= "more upreg. in late FSGS","5"= "more upreg. in early FSGS","6"= "changes in diff. directions") )
new_df$color <- revalue( as.factor(new_df$color) , c("1"= "grey","2"= "darkblue","3"= "darkgreen","4"= "magenta","5"= "orange","6"= "salmon") )

new_df <- new_df[,-c(5:8)]

new_df[["Genes"]] <- I(lapply(new_df$Row.names, function(x)  { y<- unique(entr2gname$external_gene_name[entr2gname$entrezgene%in%unlist(x)])
paste(y, collapse = ',')}))
new_df[["Genes.in.secondary"]] <- I(lapply(new_df$Row.names, function(x) { y <- unique(entr2gname$external_gene_name[entr2gname$entrezgene%in%unlist(x)])
paste(y, collapse = ',')}))


######### 2D plot for SPIA results
spia_results <- list.files(pattern = "w_SPIA_","/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/SPIA/",full.names = T)

# string wrapper
wrap_text <- function(string, n=40) 
{ require(stringr)
  if(nchar(string) > n){
    spaces <- str_locate_all(string, " ")[[1]][,1]
    chars  <- nchar(string)
    for(i in 1:floor(chars/n)) {
      s <- spaces[which.min(abs(spaces - n*i))]
      substring(string, s, s) <- "\n "
    }
    return(string)
  } else string <- string
  
}

# ECM-receptor interaction
# Cell Cycle 
# Inflammatory mediator regulation of TRP channels

# paths prepared with graphit package
path_REACT
path_KEGG

WT1_4w$entrID <- tx2gene$entrezgene[match(rownames(WT1_4w),tx2gene$ensembl_gene_id)]
WT1_4w$uniprotID <- tx2gene$uniprot_gn[match(rownames(WT1_4w),tx2gene$ensembl_gene_id)]
WT1_4w$gName <- tx2gene$external_gene_name[match(rownames(WT1_4w),tx2gene$ensembl_gene_id)]

WT1_12w$entrID <- tx2gene$entrezgene[match(rownames(WT1_12w),tx2gene$ensembl_gene_id)]
WT1_12w$uniprotID <- tx2gene$uniprot_gn[match(rownames(WT1_12w),tx2gene$ensembl_gene_id)]
WT1_12w$gName <- tx2gene$external_gene_name[match(rownames(WT1_12w),tx2gene$ensembl_gene_id)]

### add LFC of pathway members 
spia_results_LFC <- vector("list",length = length(spia_results) ) # initialize a list for results of LFC computing
# do LFC computation
{
  # if (all(spia_path$Name%in%names(path_KEGG))) {  # figure out if it's KEGG or reactome results 
  #spia_path <- spia_path[which(spia_path$Name%in%names(path_KEGG)),]
  
  kegg_path_LFC <- t(sapply( seq(path_KEGG), function(x) { pathGenes_entrID <- sub("ENTREZID:","",nodes(path_KEGG[[x]]))
  meanl2FC <-  c(mean(WT1_4w$log2FoldChange[which(WT1_4w$entrID%in%pathGenes_entrID)],na.rm=TRUE) ,mean(WT1_12w$log2FoldChange[which(WT1_12w$entrID%in%pathGenes_entrID)],na.rm=TRUE),length(nodes(path_KEGG[[x]])))
  return(meanl2FC)}))
  kegg_path_LFC <- data.frame(Name=names(path_KEGG),size=kegg_path_LFC[,3],meanl2FC_w4=kegg_path_LFC[,1],meanl2FC_w12=kegg_path_LFC[,2])
  
  reactome_path_LFC <- t(sapply( seq(path_REACT), function(x) { pathGenes_entrID <- sub("UNIPROT:","",nodes(path_REACT[[x]]))
  meanl2FC <-  c(mean(WT1_4w$log2FoldChange[which(WT1_4w$uniprotID%in%pathGenes_entrID)],na.rm=TRUE) ,mean(WT1_12w$log2FoldChange[which(WT1_12w$uniprotID%in%pathGenes_entrID)],na.rm=TRUE),length(nodes(path_REACT[[x]])))
  return(meanl2FC)}))
  reactome_path_LFC <- data.frame(Name=names(path_REACT),size=reactome_path_LFC[,3],meanl2FC_w4=reactome_path_LFC[,1],meanl2FC_w12=reactome_path_LFC[,2])
  
}

dataLFC <- kegg_path_LFC

GOplot_2Dlfc <- function(dataLFC,type="KEGG"){
  require(ggrepel)
  require(ggplot2)
  
  set1 <- readRDS(spia_results[3])
  set2 <- readRDS(spia_results[1])
  
  # union of significant pathways
  pathUnion <- union(set1$Name[which(set1$pGFWER<0.1)],set2$Name[which(set2$pGFWER<0.1)])
  
  # combine results from different weeks in one df
  XX_lfc <- data[which(!is.na(dataLFC$meanl2FC_w12)),]
  XX_lfc_SIG <- XX_lfc[which(XX_lfc$Name%in%pathUnion),]
  XX_lfc_SIG$tA_4w <- set1$tA[match(XX_lfc_SIG$Name,set1$Name)]
  XX_lfc_SIG$tA_12w <- set2$tA[match(XX_lfc_SIG$Name,set2$Name)]
  
  # Create prediction interval data frame with upper and lower lines corresponding to sequence covering minimum and maximum of x values in original dataset
  newx <- lm(as.numeric(XX_lfc_SIG$tA_12w)~as.numeric(XX_lfc_SIG$tA_4w))
  pred_interval <- predict(newx, interval="prediction", level = 0.70)
  pred_interval <- as.data.frame(pred_interval)
  
  # which labels to put
  XX_lfc_SIG$diffGO <- ifelse(XX_lfc_SIG$tA_12w < pred_interval$lwr |  XX_lfc_SIG$tA_12w > pred_interval$upr, as.character(XX_lfc_SIG$Name), "")
  # # what colors to labeled terms to use
  XX_lfc_SIG$diffGOcol <- as.factor(ifelse(XX_lfc_SIG$tA_12w < pred_interval$lwr |  XX_lfc_SIG$tA_12w > pred_interval$upr, "red", "grey"))
  
  # #  assign differnt colours to the terms based on type of annotation
  # library(RColorBrewer)
  # myColors <- c("grey",brewer.pal(4,"Set1"))
  # names(myColors) <- levels(XX_lfc$diffGOcol)
  # colScale <- scale_colour_manual(name = "grp",values = myColors) 
  
    XX_lfc_SIG$textBreaks <- sapply(as.character(XX_lfc_SIG$Name), wrap_text)
  
  ggplot(as.data.frame(XX_lfc_SIG), aes(x=as.numeric(XX_lfc_SIG$tA_4w),y=as.numeric(XX_lfc_SIG$tA_12w))) +  
    #geom_line(aes(x=XX_lfc_SIG$meanl2FC_w4,y=pred_interval$upr),linetype = 2) +
    #geom_line(aes(x=XX_lfc_SIG$meanl2FC_w4,y=pred_interval$lwr),linetype = 2) +
    geom_point(aes(size = size), colour="red",shape=21, stroke = 2) + #colScale + 
    #geom_smooth(method='lm',se =F) + # add linear regression line
    #geom_line(data=pred_interval, aes(x = fit, y = lwr), colour = "blue")+
    #geom_ribbon(data=pred_interval, aes(ymin = lwr, ymax = upr), fill = "blue", alpha = 0.2) + # add CI ribbon
    scale_size(range = c(2, 30)) +  # adjust range of bubbles sizes
    geom_text_repel(aes(label = textBreaks ),box.padding   = 0.35, point.padding = 0,segment.color = 'grey50') + #ensure that labels are not jammed together
    xlab("tA week 4") + ylab("tA week 12") + ggtitle("activity changes in KEGG pathways identified with SPIA") + 
    #scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
    theme(plot.title = element_text(face="bold",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"), 
          legend.text = element_text(size = 20),legend.title = element_text(size=20),legend.key.size = unit(2, "cm"),strip.text.y = element_text(size = 20))+
    guides(colour = guide_legend(override.aes = list(size=20))) # control size of of points in legend 
  
}



####======== Circular chord plot ==========
# description on http://wencke.github.io/ 

# # choose pathways to plot
# path_sel <- read.table(sep = "\t","/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/SPIA/SPIAKEGGREACTOMEOVERLAP_2018.csv")
# path_selKEGG <- as.character(path_sel[which(path_sel[,2]=="KEGG"),1])
# path_selREACT <- as.character(path_sel[which(path_sel[,2]=="REACTOME"),1])


### function to make a chord plot given results of SPIA and DEseq2
chord_plotLFC <- function( expr1=WT1_4w , expr2=WT1_12w ,  
                           path_selKEGG = path_selKEGG,  sortBy = 2,  
                           path_selREACT=path_selREACT, DEseq2qvalthresh=0.01 )
  {
  
  require(graphite)
  require("qdapTools")
  
  # load necessary functions 
  source("/data/user/tpadvits/PROJECTS/Helper_chordPlot.R")
  
  # load graphite pathway data
  load("/data/user/tpadvits/PROJECTS/MARTIN_PJ/RNAseq_analysis/NPHS2/SPIA/SPIApathInfo.02.05.19.rda")
  
    SELreact_LFC <- lapply( path_selREACT , function(x) tryCatch( nodes(path_REACT[[x]]),error=function(e) numeric() ) )
    names(SELreact_LFC) <- path_selREACT
    SELkegg_LFC <- lapply( path_selKEGG , function(x)  tryCatch(nodes(path_KEGG[[x]]),error=function(e) numeric() ) )
    names(SELkegg_LFC) <- path_selKEGG
  

  # # choose significant patwhays
  # all_SPIA_REACT_sig <- lapply(all_SPIA_REACT, function(x) subset(x, x$pGFdr< pGFDRthresh))
  # all_SPIA_KEGG_sig <- lapply(all_SPIA_KEGG, function(x) subset(x, x$pGFdr< pGFDRthresh))
  # 
  
  ### LFC of genes from Selected KEGG and Reactome pathways 
  # convert lists to df
  SELreact_LFC_df <- list2df(SELreact_LFC)
  SELkegg_LFC_df <-  list2df(SELkegg_LFC)
  
  # convert specific IDs to gene names
  SELreact_LFC_df[,1] <- tx2gene$external_gene_name[match(sub("ENTREZID:","",SELreact_LFC_df[,1]),tx2gene$entrezgene_id)] 
  SELkegg_LFC_df[,1] <- tx2gene$external_gene_name[match(sub("ENTREZID:","",SELkegg_LFC_df[,1]),tx2gene$entrezgene_id)]
  
  # remove duplicate rows and NA genes 
  SELreact_LFC_df <- SELreact_LFC_df[!duplicated(SELreact_LFC_df),]
  SELreact_LFC_df <- SELreact_LFC_df[!is.na(SELreact_LFC_df[,1]),]
  SELkegg_LFC_df <- SELkegg_LFC_df[!duplicated(SELkegg_LFC_df),]
  SELkegg_LFC_df <- SELkegg_LFC_df[!is.na(SELkegg_LFC_df[,1]),]
  
  # add LFC
  SELreact_LFC_df$w4lfc <- expr1$log2FoldChange[ match(SELreact_LFC_df$X1,expr1$gName)]
  SELreact_LFC_df$w4qval <-  expr1$padj[ match(SELreact_LFC_df$X1,expr1$gName) ]
  SELreact_LFC_df$w12lfc <- expr2$log2FoldChange[ match(SELreact_LFC_df$X1,expr2$gName)]
  SELreact_LFC_df$w12qval <-  expr2$padj[ match(SELreact_LFC_df$X1,expr2$gName)]
  SELreact_LFC_dfSig <- subset( SELreact_LFC_df, SELreact_LFC_df$w4qval< DEseq2qvalthresh |  SELreact_LFC_df$w12qval< DEseq2qvalthresh )
  SELreact_LFC_dfSig$X2 <- paste( SELreact_LFC_dfSig$X2,"REACTOME", sep = "^")
  
  SELkegg_LFC_df$w4lfc <- expr1$log2FoldChange[match(SELkegg_LFC_df$X1,expr1$gName)]
  SELkegg_LFC_df$w4qval <-  expr1$padj[match(SELkegg_LFC_df$X1,expr1$gName)]
  SELkegg_LFC_df$w12lfc <- expr2$log2FoldChange[match(SELkegg_LFC_df$X1,expr2$gName)]
  SELkegg_LFC_df$w12qval <-  expr2$padj[match(SELkegg_LFC_df$X1,expr2$gName)]
  SELkegg_LFC_dfSig <- subset(SELkegg_LFC_df, SELkegg_LFC_df$w4qval< DEseq2qvalthresh |  SELkegg_LFC_df$w12qval< DEseq2qvalthresh )
  SELkegg_LFC_dfSig$X2 <- paste(SELkegg_LFC_dfSig$X2,"KEGG", sep = "^")
  
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
  print(pp)
  
  return( SELboth_LFC_dfSig)
}


  


####======== Heatmaps of important genes ==========
# gene targets of WT1hetdel mouse
WT1_closGene <- read.table( sep = "\t" , header = F , file="/data/user/tpadvits/PROJECTS/MARTIN_PJ/DATA_integ/WT1_Tead1/coocurance/WT1_nearestg.bed")
WT1_closGene_sig <- unique (sub("\\..*","",WT1_closGene$V10[which(abs(WT1_closGene$V13)<=2000)]))
WT1_closGene_sig_genes <- unique(tx2gene$external_gene_name[which(tx2gene$ensembl_transcript_id%in%(WT1_closGene_sig))])


### interesting SPIA pathways
spia_results <- list.files(pattern = "w_SPIA_","/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/SPIA/",full.names = T)
ephr_paths_react <- do.call( union , Filter(length,lapply(spia_results, function(x) grep("ephr",readRDS(x)$Name,value = T,ignore.case = T))))
coll_paths_react <- do.call( union , Filter(length,lapply(spia_results, function(x) grep("collag",readRDS(x)$Name,value = T,ignore.case = T))))


###============= heatmap of the expression across replicates
# read expression data
expression_dat <- read.table(header = T, stringsAsFactors = F, row.names = 1, "/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/RNAseq/DE_analysis/expr_IntExon_genes_tpmSized.csv")
colnames(expression_dat) =c("w4_het.del","w4_wt","w4_wt","w4_het.del","w4_wt","w4_het.del","w12_wt","w12_wt","w12_het.del","w12_het.del","w12_het.del","w12_wt")
expression_dat$uniprotID <- tx2gene$uniprot_gn[match(rownames(expression_dat),tx2gene$ensembl_gene_id)]
expression_dat$gName <- tx2gene$external_gene_name[match(rownames(expression_dat),tx2gene$ensembl_gene_id)]
expression_dat_filt <- expression_dat[which(expression_dat$gName %in% WT1_4w$gName),]

### Ephrin and collagene pathways
{   
  ephr_react_TPM <- lapply( ephr_paths_react, function(x) { pathGenes_entrID <- sub("UNIPROT:","",nodes(path_REACT[[x]]))
    l2FC <- expression_dat_filt[which(expression_dat_filt$uniprotID%in%pathGenes_entrID),1:12]
    rownames(l2FC) <-  expression_dat_filt$gName[which(expression_dat_filt$uniprotID%in%pathGenes_entrID)]
    return(l2FC) })
  ephr_react_TPM <- do.call(rbind, ephr_react_TPM)
  ephr_react_TPM <- ephr_react_TPM[!duplicated(ephr_react_TPM),]
  
  collag_react_TPM <- lapply( coll_paths_react, function(x) { pathGenes_entrID <- sub("UNIPROT:","",nodes(path_REACT[[x]]))
    l2FC <-  expression_dat_filt[which(expression_dat_filt$uniprotID%in%pathGenes_entrID),1:12]
    rownames(l2FC) <-  expression_dat_filt$gName[which(expression_dat_filt$uniprotID%in%pathGenes_entrID)]
    return(l2FC) })
  collag_react_TPM <- do.call(rbind, collag_react_TPM)
  collag_react_TPM <- collag_react_TPM[!duplicated(collag_react_TPM),]
  
} 
### collagene genes 
{ collagene_TPM <- read.table(sep = "\t", stringsAsFactors = F, row.names = 1, header = T, "/data/user/tpadvits/PROJECTS/MARTIN_PJ/Publication/WT1_paper/expr_heatmaps/collagens_expr/collagen_tpm.csv")
  collagene_TPM <- collagene_TPM[which(rownames(collagene_TPM)%in%expression_dat_filt$gName),]
  collagene_TPM_ave <- cbind(rowMeans(collagene_TPM[,1:3]),rowMeans(collagene_TPM[,4:6]),rowMeans(collagene_TPM[,7:9]),rowMeans(collagene_TPM[,10:12]))
  colnames(collagene_TPM_ave) <- c("w4_wt","w4_het.del","w12_wt","w12_het.del")
}

library(ggplot2)
library(reshape)
#XX <- t(scale(t(ephr_react_LFC),center = F))
XX <- collagene_TPM_ave[which(rowMeans(collagene_TPM_ave) > 1),]
XX <- t(scale(t(XX[,1:4]),center = F))
ord <- hclust( dist(XX, method = "euclidean"), method = "ward.D" )$order
XX.m <- melt( as.matrix( XX ) )
XX.m$X1 <- factor( XX.m$X1, levels = rownames(XX)[ord], labels =  rownames(XX)[ord] )


ggplot(XX.m, aes(X2, X1)) + 
  geom_tile(aes(fill = value), colour = "darkgrey") + theme_bw() + theme( panel.grid.major = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
                                                                         panel.border = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
                                                                                             axis.text.x = element_text(face="bold", angle=45,vjust=0.5,size=12),axis.text.y = element_text(size=10,face="italic",family = "arial"))+ 
  # geom_text(aes(label = round(padjust, 4))) +    axis.text.x = element_text(face="bold", angle=0,vjust=0.5,size=11),axis.text.y = element_text(size=11,face="italic",family = "arial"))
  scale_fill_gradient2( midpoint = 0.8, low="white" ,mid = "yellow", high = "red", na.value = "grey",name="scaled tpm of gene expression") 
  # scale_fill_manual(values=my_palette, breaks=c(levels(as.factor(XX$value))[seq(1, 160, by=10)],2.979),name="log2 Fold.Change",na.value = "red")

###================ LFC heatmap 
### Ephrin and collagene pathways
{ 
  ephr_react_LFC <- lapply( ephr_paths_react, function(x) { pathGenes_entrID <- sub("UNIPROT:","",nodes(path_REACT[[x]]))
  l2FC <-  cbind( WT1_4w$log2FoldChange[which(WT1_4w$uniprotID%in%pathGenes_entrID)] , WT1_12w$log2FoldChange[which(WT1_12w$uniprotID%in%pathGenes_entrID)] )
  rownames(l2FC) <-  WT1_4w$gName[which(WT1_4w$uniprotID%in%pathGenes_entrID)]
  return(l2FC) })
  ephr_react_LFC <- do.call(rbind, ephr_react_LFC)
  ephr_react_LFC <- ephr_react_LFC[!duplicated(ephr_react_LFC),]
  
  collag_react_LFC <- lapply( coll_paths_react, function(x) { pathGenes_entrID <- sub("UNIPROT:","",nodes(path_REACT[[x]]))
  l2FC <-  cbind( WT1_4w$log2FoldChange[which(WT1_4w$uniprotID%in%pathGenes_entrID)] , WT1_12w$log2FoldChange[which(WT1_12w$uniprotID%in%pathGenes_entrID)] )
  rownames(l2FC) <-  WT1_4w$gName[which(WT1_4w$uniprotID%in%pathGenes_entrID)]
  return(l2FC) })
  
  collag_react_LFC <- do.call(rbind, collag_react_LFC)
  collag_react_LFC <- collag_react_LFC[!duplicated(collag_react_LFC),]
  
  colnames(ephr_react_LFC)   <- colnames(collag_react_LFC) <- c("week 4", "week 12")
  }

### collagene genes 
collagene_l2FC <-  cbind( WT1_4w[which(WT1_4w$gName%in%rownames(collagene_TPM)),c("log2FoldChange","gName")] , WT1_12w[which(WT1_12w$gName%in%rownames(collagene_TPM)),c("log2FoldChange","gName")] )
rownames(collagene_l2FC) <- collagene_l2FC$gName
collagene_l2FC <- collagene_l2FC[,c(1,3)]
colnames(collagene_l2FC) <- c("week4","week12")


library(ggplot2)
library(reshape)
#XX <- t(scale(t(ephr_react_LFC),center = F))
XX <- collagene_l2FC
#XX <- t(scale(t(collagene_l2FC),center = F))
ord <- hclust( dist(XX, method = "euclidean"), method = "ward.D" )$order
XX.m <- melt( as.matrix( XX ) )
XX.m$X1 <- factor( XX.m$X1, levels = rownames(XX)[ord], labels =  rownames(XX)[ord] )


ggplot(XX.m, aes(X2, X1)) + 
  geom_tile(aes(fill = value), colour = "darkgrey") + theme_bw() + theme( panel.grid.major = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
                                                                          panel.border = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
                                                                          axis.text.x = element_text(face="bold", angle=45,vjust=0.5,size=12),axis.text.y = element_text(size=10,face="italic",family = "arial"))+ 
  # geom_text(aes(label = round(padjust, 4))) +    axis.text.x = element_text(face="bold", angle=0,vjust=0.5,size=11),axis.text.y = element_text(size=11,face="italic",family = "arial"))
  scale_fill_gradient2( low="blue" ,mid = "white", high = "red", na.value = "grey",name="scaled log2 fold change\nof gene expression") 
# scale_fill_manual(values=my_palette, breaks=c(levels(as.factor(XX$value))[seq(1, 160, by=10)],2.979),name="log2 Fold.Change",na.value = "red")

 