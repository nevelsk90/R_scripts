#### FUNCTIONS ####

### define function for making 1) sample distance heatmap and 2) PCA 
PCA_heatm_plot=function(count_table, sampleNames=NA,
                        groupsExp , PCs = c(1,2),
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
  require(matrixStats)
  
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
  print( summary(pr_comp_rnorm))
  #biplot(pr_comp_rnorm)
  XX = as.data.frame(pr_comp_rnorm$x)
  g <- qplot( x=XX[,PCs[1]], y=XX[,PCs[2]], 
              data=XX,  colour=factor( groupsExp ), main=title) + geom_point(size=4 ) +
    xlab(paste("PC", PCs[1], sep=""))+ylab(paste("PC", PCs[2], sep=""))
  # labels and dots sizes
  g+theme_bw() + theme(text=element_text(size=20))+
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
fun_homoTO.FROMmouse <- function(gns, TO=T ){
  # if TO parameter is TRUE then Homo -> Mouse, if FALSE then Mouse -> Homo
  
  options(connectionObserver = NULL)
  require( "org.Hs.eg.db" )
  require( "org.Mm.eg.db")
  require( "Orthology.eg.db")
  require("AnnotationDbi")
 
  if(isTRUE(TO)){
    egs <- mapIds(org.Hs.eg.db, gns, "ENTREZID","SYMBOL")
    mapped <- AnnotationDbi::select(Orthology.eg.db, egs, "Mus.musculus","Homo.sapiens")
    mapped$MUS <- mapIds(org.Mm.eg.db, as.character(mapped$Mus.musculus), "SYMBOL", "ENTREZID")
    mapped$MUS.ens <- mapIds(org.Mm.eg.db, as.character(mapped$Mus.musculus), "ENSEMBL", "ENTREZID")
  } else {
    egs <- mapIds(org.Mm.eg.db, gns, "ENTREZID","SYMBOL")
    mapped <- select( Orthology.eg.db, egs, "Homo.sapiens", "Mus.musculus")
    mapped$HOMO <- mapIds(org.Hs.eg.db, as.character(mapped$Homo.sapiens), "SYMBOL", "ENTREZID")
    mapped$HOMO[ is.na(mapped$Homo.sapiens)] <-  toupper(rownames(mapped))[ is.na( mapped$Homo.sapiens ) ] 
    mapped$HOMO[rownames(mapped)=="Cd59a"]<- "CD59"
    mapped$HOMO.ens <- mapIds(org.Hs.eg.db, as.character(mapped$HOMO), "ENSEMBL", "SYMBOL")
    
  }
  
  return(mapped)
}


### function to sum P-values using Fisher method
fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)


### moving average  func
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

### wrap text strins
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


#### https://stackoverflow.com/questions/15006211/how-do-i-generate-a-mapping-from-numbers-to-colors-in-r
### numeric to color
map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

### Range standardization (0 to 1) in R  https://stackoverflow.com/questions/5665599/range-standardization-0-to-1-in-r 
range01 <- function(x){(x-min(x))/(max(x)-min(x))}


#### functional annotation by clusterprofile

cProfiler.GKR <-  function( ggenes.eID, 
                            gene.bckgrnd_eID ){
  ## ggenes.eID and gene.bckgrnd_eID are entrez IDs 
  
  require(clusterProfiler)
  require(org.Mm.eg.db)
  require(ReactomePA)
  
  GO.enrich <- clusterProfiler::enrichGO( gene=ggenes.eID , 
                                          OrgDb = "org.Mm.eg.db",
                                          ont="ALL",
                                          universe= gene.bckgrnd_eID,
                                          minGSSize =3 ,
                                          readable = T
  )
  KEGG.enrich <- clusterProfiler::enrichKEGG( gene= ggenes.eID , 
                                              organism = "mmu",
                                              universe= gene.bckgrnd_eID, 
                                              minGSSize =3 ,
                                              use_internal_data = T
  )
  KEGG.enrich@result <- cbind( ONTOLOGY= "KEGG" , KEGG.enrich@result  ) 
  
  REACT.enrich<- ReactomePA::enrichPathway(gene=ggenes.eID , 
                                           organism = "mouse",
                                           universe = gene.bckgrnd_eID,
                                           minGSSize =3 ,
                                           readable=T )
  REACT.enrich@result <- cbind( ONTOLOGY= "REACT" , REACT.enrich@result  ) 
  
  ll <- list(GO.enrich, KEGG.enrich ,REACT.enrich)          
  names(ll)<- c("GO.enrich", "KEGG.enrich" ,"REACT.enrich")  
  return(ll)
}


### calculate AUCell score for a group of gene sets
## the function can handle both named numeric and charachter vectors

calculateAUCcell <- function( geneSets , AUCthrsh=0.05, ddata ,
                              Seurat.type =F )
  {
  require(GSEABase)
  ## the function can handle both named numeric and charachter vectors
  
  # split genesets by direction of change if numeric vectors are provided
  if( is.numeric(geneSets[[1]])) {
    geneSets.Split <- Reduce( c , lapply(seq( geneSets), function(jj){
      UP <- names(geneSets[[jj]])[geneSets[[jj]]>0]
      DOWN <- names(geneSets[[jj]])[geneSets[[jj]]<0]
      return(list(UP,DOWN))
    }))
    names( geneSets.Split ) <- c(rbind(  paste( names(geneSets) ,"UP",sep = "_"),
                                         paste( names(geneSets) ,"DOWN",sep = "_") ))

  } else geneSets.Split <- geneSets
  
  # prepare gene sets
  genesets_list <- lapply( seq(geneSets.Split) , function( ii ) {
    GeneSet( geneSets.Split[[ii]], 
             setName= names(geneSets.Split)[ii] ) } )
  genesets_list <-  GSEABase::GeneSetCollection( genesets_list )
  
  # load expression data 
  if( isTRUE(Seurat.type) ) {
    exprMatrices <- ddata@assays$RNA@counts
  }  else   exprMatrices <- as.matrix(ddata)
 

  # 1. Build gene-expression rankings for each cell  
  cells_rankings <-  AUCell::AUCell_buildRankings( exprMatrices, 
                                                   nCores=1, plotStats=F)
  
  # 2. Calculate enrichment for the gene signatures (AUC)
  # store NAs if less than 20% of geneSet genes are expressed in snRNAseq
  cells_AUC <- AUCell::AUCell_calcAUC( genesets_list , cells_rankings ,
                                       aucMaxRank= ceiling(AUCthrsh * nrow(cells_rankings)), 
                                       verbose = T )
  cells_AUC <- AUCell::getAUC( cells_AUC)
  
  # subtract down scores from the UP scores
  if( is.numeric(geneSets[[1]])) { 
    cells_AUCres <- t(cells_AUC)
    cells_AUCres <- cells_AUCres[,seq(1,ncol(cells_AUCres),2)] - cells_AUCres[,seq(2,ncol(cells_AUCres),2)]
    colnames(cells_AUCres) <- names(geneSets)
    cells_AUCres <- t(cells_AUCres)
    } else cells_AUCres <- cells_AUC
  
  return(cells_AUCres)
}

### https://stackoverflow.com/questions/37613345/r-convert-upper-triangular-part-of-a-matrix-to-symmetric-matrix
ultosymmetric_diagonalone=function(m){
  m = m + t(m) - 2*diag(diag(m)) + diag(1,nrow=dim(m)[1])
  return (m)}

makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

## https://stackoverflow.com/questions/13112238/a-matrix-version-of-cor-test
cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}


###https://stats.stackexchange.com/questions/160671/estimate-lag-for-granger-causality-test
select.lags<-function(x,y,max.lag=8) {
  y<-as.numeric(y)
  y.lag<-embed(y,max.lag+1)[,-1,drop=FALSE]
  x.lag<-embed(x,max.lag+1)[,-1,drop=FALSE]
  
  t<-tail(seq_along(y),nrow(y.lag))
  
  ms=lapply(1:max.lag,function(i) lm(y[t]~y.lag[,1:i]+x.lag[,1:i]))
  
  pvals<-mapply(function(i) anova(ms[[i]],ms[[i-1]])[2,"Pr(>F)"],max.lag:2)
  ind<-which(pvals<0.05)[1]
  ftest<-ifelse(is.na(ind),1,max.lag-ind+1)
  
  aic<-as.numeric(lapply(ms,AIC))
  bic<-as.numeric(lapply(ms,BIC))
  structure(list(ic=cbind(aic=aic,bic=bic),pvals=pvals,
                 selection=list(aic=which.min(aic),bic=which.min(bic),ftest=ftest)))
}


### Transformation of GRange object to density per bin
# https://divingintogeneticsandgenomics.com/post/compute-averages-sums-on-granges-or-equal-length-bins/
averagePerBin <- function(x, binsize, mcolnames=NULL)
{
  if (!is(x, "GenomicRanges"))
    stop("'x' must be a GenomicRanges object")
  if (any(is.na(seqlengths(x))))
    stop("'seqlengths(x)' contains NAs")
  bins <- IRangesList(lapply(seqlengths(x),
                             function(seqlen)
                               IRanges(breakInChunks(seqlen, binsize))))
  ans <- as(bins, "GRanges")
  seqinfo(ans) <- seqinfo(x)
  if (is.null(mcolnames))
    return(ans)
  averageMCol <- function(colname)
  {
    cvg <- coverage(x, weight=colname)
    views_list <- RleViewsList(
      lapply(names(cvg),
             function(seqname)
               Views(cvg[[seqname]], bins[[seqname]])))
    unlist(viewMeans(views_list), use.names=FALSE)
  }
  mcols(ans) <- DataFrame(lapply(mcols(x)[mcolnames], averageMCol))
  ans
}

### How to prevent reduce when using setdiff on GRanges  https://www.biostars.org/p/489350/ 
# Written by bosberg on Biostars, March 22, 2021. Free to use for scholarly purposes.
GRanges_subtract <- function( gr1, gr2 )
{
  require(bedtoolsr)
  require(GenomicRanges)
  # Subtract GRange object gr2 from gr1, but unlike setdiff, preserve individual
  # ranges in gr1
  df_1 = data.frame( seqnames=seqnames(gr1), start=start(gr1)-1, end=end(gr1), strand=strand(gr1), mcols( gr1 ) )
  df_2 = data.frame( seqnames=seqnames(gr2), start=start(gr2)-1, end=end(gr2), strand=strand(gr2), mcols( gr2 ) )
  #                                                          ^ -1 --> convert to base-0 start for bedtools
  result = bedtoolsr::bt.subtract(df_1, df_2)
  
  if ( length(result)==0 ){
    # subtraction has left nothing remaining. Return empty GRanges obj.
    return( GRanges() )
  } else {
    
    colnames( result ) = colnames( df_1 )
    result$start=result$start+1
    #                        ^ reset to base-1 notation consistent with GRanges
    
    return ( GRanges( result ) )
  }
}


#random drawings from a custom probability distribution derived from a dataset
# https://stackoverflow.com/questions/51500331/r-random-drawings-from-a-custom-probability-distribution-derived-from-a-datase 
mysampler <- function(x, n) {
  y <- runif(n)
  c(quantile(x, probs = y)) #maybe a different quantile type should be used?
}

## compute the size of directory in R
## https://stackoverflow.com/questions/39753172/compute-the-size-of-directory-in-r
dir_size <- function(path, recursive = TRUE) {
  stopifnot(is.character(path))
  files <- list.files(path, full.names = T, recursive = recursive)
  vect_size <- sapply(files, function(x) file.size(x))
  size_files <- sum(vect_size)
  size_files
}

### Convert read counts to transcripts per million (TPM).
# https://gist.github.com/slowkow/c6ab0348747f86e2748b 
counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

###Jaccard similarity
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
