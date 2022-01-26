#### FUNCTIONS ####

### define function for making 1) sample distance heatmap and 2) PCA 
PCA_heatm_plot=function(count_table, sampleNames=NA,
                        groupsExp ,
                        topNprct=0.1 ,
                        logtrans=F, title="PCA")
{
  # sampleNames - a charachter vector to use as labels of individual samples
  # groupsExp - a charachter vector for experimental groups, e.g. wild-type and mutant
  # logtrans - logic value, whether log transform the data or not 
  # topNprct - proportions of genes with highest variance to consider for PCA, 0.1 by default
  
  ### Estimate and plot eucledean disances between samples
  # !!!! remember that not TPM but count data shall be normalized with rlog !!!!
  library("RColorBrewer")
  library("pheatmap")  
  require(ggplot2)
  
  # define sample names as column names if sample names are not provided
  if ( length(sampleNames) == ncol(count_table)) colnames(count_table) <- sampleNames else sampleNames <- colnames(count_table)
  
  # select top {topNprct} portion of most variable genes
  count_table <- as.matrix(count_table)
  count_table <- count_table[ which( rowVars(count_table) > quantile( rowVars(count_table) , 1-topNprct )) ,  ]
  
  # log transform if necessary
  if( isTRUE(logtrans)) rld <- log( count_table + 1 ) else rld <- count_table
  rownames(rld)=rownames(count_table)
  
  # calculate distances
  euclDists <- dist( t( rld ) )
  euclDistsMatrix <- as.matrix( euclDists )
  rownames(euclDistsMatrix) <- names(count_table)
  colnames(euclDistsMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(euclDistsMatrix,
           clustering_distance_rows=euclDists,
           clustering_distance_cols=euclDists,
           col=colors , labels_row =  sampleNames )
  
  
  ## PCA for count data rlog transform  using DEseq2 algorythm
  pr_comp_rnorm <- prcomp(t(rld),scale=F)
  #biplot(pr_comp_rnorm)
  XX = as.data.frame(pr_comp_rnorm$x)
  g <- qplot( x=PC1, y=PC2, data=XX,  colour=factor(groups), main=title) + 
    geom_point(size=4  ) 
  # labels and dots sizes
  g+theme_bw() + theme(axis.text=element_text(size=20),
                       axis.title=element_text(size=24,face="bold"), legend.text = element_text(size = 24),legend.title = element_text(size=24))+
    scale_size(range = c(5, 6))+theme(legend.key.size = unit(2, "cm"))+ labs(colour = "condition")+
    geom_text(aes(label=rownames(XX)),hjust=0.4, vjust=-0.8,size=5)
}

# ## file for UCSC <> ENSEMBLE chromosome name conversion
# chrUCSCtoENSEMBLE <- read.table("/data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38/ensembTOucsc_chromname.txt",sep="\t")

### TFtarget caler from bed file
TFtgClosGene_func <- function( bedFile , glvl=F ,
                               Gene_annot_forTFcaller = "")
  {
  require(TFTargetCaller)
  peaks <- read.table( file=bedFile , sep = "\t", header = T)
  peaks <- dplyr::mutate(peaks,center= round(rowMeans(peaks[,c("start","stop")])))
  peaks=peaks[, c( 1, ncol(peaks) )]
  names(peaks)=c("chromosome", "center")
  peaks$chromosome=as.character(peaks$chromosome)
  

  ## do tgenes inference
  tgene_closGene_qvalue <- TFTargetCaller( peaks, Gene_annot_forTFcaller, method="ClosestGene", ClosestGeneScore="qvalue")
  tgene_closGene_score <- TFTargetCaller( peaks, Gene_annot_forTFcaller, method="ClosestGene")
  tgene_closGene=as.data.frame(cbind(tgene_closGene_qvalue ,tgene_closGene_score))
  
  if (glvl==T){
    # summarize to gene lvl only significant
    ClosGene_glvl <- tgene_closGene[which(rownames(tgene_closGene)%in%tx2gene$ensembl_transcript_id & tgene_closGene[,1]<0.1),]
    ClosGene_glvl <- merge(ClosGene_glvl,tx2gene[,1:2],by.x = 0 , by.y = "ensembl_transcript_id")
    colnames(ClosGene_glvl) <- c("trscrID","qvalue","score","geneID")
    # aggregate by score
    ClosGene_glvl <- ClosGene_glvl[order(ClosGene_glvl$geneID, -abs(ClosGene_glvl$score) ),]
    tgene_closGene <- aggregate(score~geneID,data=ClosGene_glvl,FUN=sum)
  }
  
  
  return(tgene_closGene)
}

### Function to import and summarize transcript-level abundance estimates 
### from kallisto output to gene-level estimates
### and then summarise Exonic and Intronic results
kallistoExInt_sum <- function( dir_exo , dir_intro , samples_exo , samples_intro, 
                               nameE, tx2gene=tx2gene ) 
  {
  require("tximport")
  require("biomaRt")
  require(readr)
  require(rhdf5)
  
  ### get a table with different IDs and annotation from the BioMart
  tx2gene <- dplyr::rename( tx2gene, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id )
  
  ### if transcript IDs from kallisto output has decimal (mouse)
    if (length(grep("*\\..*",
                    read.table( file=kall_files_exons[1],
                    sep = "\t",header=T)[,1])) !=0 ) stop("check transcript IDs in the kallisto output, they probably have decimals.\n
                                                          Change tx2gene object so that transcript id includes version.")
  
  ### EXONS
  kall_files_exons <- list.files(dir_exo,pattern = "\\.tsv",recursive = TRUE,full.names = T)
  ## sum to gene level
  txi <- tximport( kall_files_exons, type = "kallisto", tx2gene = tx2gene, importer = read_tsv ,  ignoreTxVersion = T )
  colnames(txi$counts)=samples_exo
  colnames(txi$abundance)=samples_exo
  expr_table=txi$counts
  estim_exon=txi$abundance
  
  
  ## filter non-expressed genes
  expr_table_genes_non0=expr_table[rowSums(expr_table)>=ncol(expr_table),]
  expr_table_genes_non0=round(expr_table_genes_non0)
  
  ### INTRONS
  kall_files_intro <- list.files( dir_intro,pattern = "\\.tsv",recursive = TRUE,full.names = T)
  ## sum to gene level
  txi <- tximport( kall_files_intro , type = "kallisto", tx2gene = tx2gene, importer = read_tsv  )
  colnames(txi$counts)=samples_intro
  colnames(txi$abundance)=samples_intro
  expr_intron=txi$counts
  estim_intron=txi$abundance
  
  
  ### SUM OF EXON & INTRON
  x1=cbind.data.frame(expr_table, gene_id=rownames(expr_table))
  x2=cbind.data.frame(expr_intron, gene_id=rownames(expr_intron))
  expr_IntExon=aggregate(.~gene_id, rbind(x1,setNames(x2,names(x1))), sum)
  rownames(expr_IntExon)=expr_IntExon$gene_id
  expr_IntExon=expr_IntExon[,-1]
  
  saveRDS( estim_exon, file=paste( outdir , nameE , "_Exon_TPMglvl.rda" , sep = ""))
  saveRDS( estim_intron, file=paste( outdir , nameE , "_Int_TPMglvl.rda", sep = "") )
  saveRDS( expr_IntExon, file=paste( outdir , nameE , "_IntExon_glvl.rda", sep = "") )
                             
}

## summarise to glvl seperately
trlvlTOglvl_func <- function( tgene_closGene , gName=T , thrshld=0.1)
  {
  # summarize to gene lvl only significant
  
  ClosGene_glvl <- tgene_closGene[which(rownames(tgene_closGene)%in%tx2gene$ensembl_transcript_id & tgene_closGene[,1]<thrshld),]
  
  if (gName==T){
    ClosGene_glvl <- merge(ClosGene_glvl,tx2gene[,c("ensembl_transcript_id","external_gene_name")],by.x = 0 , by.y = "ensembl_transcript_id")
    
  } else     ClosGene_glvl <- merge(ClosGene_glvl,tx2gene[,c("ensembl_transcript_id","ensembl_gene_id")],by.x = 0 , by.y = "ensembl_transcript_id")
  
  
  colnames(ClosGene_glvl) <- c("trscrID","qvalue","score","geneID")
  # aggregate by score
  ClosGene_glvl <- ClosGene_glvl[order(ClosGene_glvl$geneID, -abs(ClosGene_glvl$score) ),]
  ClosGene_glvl <- aggregate(score~geneID,data=ClosGene_glvl,FUN=sum)
  
  return(ClosGene_glvl)
}

### convert Human to mouse genes
convertHumanGeneList <- function(x ,  uniqueG =T)
  {
  
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  require("biomaRt")
  
  genesV2 = getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = x , mart = human, attributesL = c("external_gene_name"), martL = mouse, uniqueRows= uniqueG)
  
  # Print the first 6 genes found to the screen
  print(head(genesV2))
  return(genesV2)
}

### function to sum P-values using Fisher method
fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)


### moving average func
slideFunct <- function(data, window, step, Fmedian=F)
  {
  total <- length(data)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    if( isTRUE(Fmedian) ){
      result[i] <- median(data[spots[i]:(spots[i]+window)], na.rm = T)
    } else result[i] <- mean(data[spots[i]:(spots[i]+window)], na.rm = T)
  }
  return(result)
}
### compute the disease score using marker genes and scRNAseq
# disease_score <- function( scRNAseq = scPodo_libNorm , wt_cells=c(1:587) ,marker_gene_list , marker_gene_LFC , mainPlot="disease score based on scRNAseq 385DE genes")
#   {
#   scBerlin_PodoDE <- scRNAseq[ match( marker_gene_list , rownames(scRNAseq) ) , ]
#   
#   
#   PodoGene_WTmean <- rowMeans( scBerlin_PodoDE[,wt_cells] )
#   
#   PodoGene_disStat <- sapply( 1:nrow(scBerlin_PodoDE) , function(x) 
#     if ( marker_gene_LFC[x] > 0)  scBerlin_PodoDE[x,] - PodoGene_WTmean[x] else (scBerlin_PodoDE[x,] - PodoGene_WTmean[x])*( -1 ) ) 
#   
#   PodoGene_disVect <- apply( PodoGene_disStat , 1, function(x) mean(as.numeric(x) , na.rm = T) )
#   
#   # # plot disease status for every cell
#   # plot(PodoGene_disVect, xlab ="cells" , ylab="disease score", main=mainPlot , col = c(rep("blue", 587) , rep("red" , 273)))
#   # legend("top", legend=c("wild-type", "WT1 het.del."), fill=c("blue", "red"),  title = "genotype")
#   
#   return(PodoGene_disVect)
#   
# }

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
