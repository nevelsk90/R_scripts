#### calculate TFA using generalised inversed of the prior ####
# exprMat -  expression matrix of normalised sc/snRNAseq data: 
# rows are genes cells are columns
# ATACseq_prior - GRN matrix where columns are TFs which activities we want to estimate
# rows are possible targets of these TFs
TFAginv_50ss_func <- function( exprMat  ,
                               ATACseq_prior ,
                               weighted=T )
  {
    require(abind)
    
    # from the prior select only TFs expressed in podocytes and genes expressed in snRNAseq
    ATACseq_prior_expr <- ATACseq_prior[ rownames(ATACseq_prior) %in% rownames(exprMat), ]
    ATACseq_prior_expr <- ATACseq_prior_expr[rowSums(ATACseq_prior_expr)>0,]
    
    
    ### for the sake of stability run TFA estimation on differnt subsets of genes
    TFA_ginv_50ss <- lapply( 1:50, 
                             function(i, prior=ATACseq_prior_expr )
                             {
                               # choose a subset of genes that will participate in TFA estimation in this round
                               cat( i, "\n")
                               Idv <- sample(1:nrow(prior),  nrow(prior)*0.63)
                               prior_ss50 <- prior[Idv,]
                               
                               if( isFALSE(weighted)) prior_ss50[prior_ss50!=0] <- 1 
                               # targets <- names(prior_ss50[prior_ss50[,j]>0,j])
                               # to  estimate TFA use generalized inverse
                               TFA_ginv <-MASS::ginv(prior_ss50) %*% as.matrix( exprMat[match(rownames(prior_ss50),
                                                                                              rownames(exprMat)),] )
                               rownames(TFA_ginv) <- (colnames(prior_ss50))
                               
                               print(summary(TFA_ginv[,1]))
                               
                               return(TFA_ginv)
                               gc()
                             })
    
    a <- do.call(abind, c(TFA_ginv_50ss, list(along=3)))
    TFA_ginv_50ss <- apply(a, 1:2, median, na.rm = T )
    
    # rownames(TFA_ginv_50ss) <- colnames(ATACseq_prior_expr)
    # colnames(TFA_ginv_50ss) <- colnames(exprMat_list)
    
    return(TFA_ginv_50ss)
    gc()
  } 


#### calculate TFA by using average expression of targets ####
TFAave_50ss_func <- function( exprMat,
                              ATACseq_prior ,
                              weighted = T )
  {
  ## weighted argument controls if the prior matrix is binary or have weights
    require(abind)
    
    # from the prior select only TFs expressed in podocytes and genes expressed in snRNAseq
    ATACseq_prior_expr <- ATACseq_prior[ rownames(ATACseq_prior) %in% rownames(exprMat), ]
    ATACseq_prior_expr <- ATACseq_prior_expr[rowSums(ATACseq_prior_expr)>0,]
    
    ### for the sake of stability run TFA estimation on differnt subsets of genes
    TFA_average_50ss <- lapply( 1:50, function(i, prior=ATACseq_prior_expr )
      {
      # choose a subset of genes that will participate in TFA estimation in this round
      cat( i, "\n")
      Idv <- sample(1:nrow(prior),  nrow(prior)*0.63)
      prior_ss50 <- prior[Idv,]
      
      # to  estimate TFA simply use average of TF targets
      TFA_mean <- Reduce( rbind, 
                          lapply( 1:ncol(prior) , 
                                  function(j)
                                  {
                                    targets <- rownames(prior_ss50)[prior_ss50[,j]>0]
                                    if( isTRUE( weighted )){ 
                                      out <- tryCatch( colMeans( prior_ss50[,j] * exprMat[ match( rownames(prior_ss50) , rownames(exprMat) ) , ] , na.rm = T) ,
                                                       error=function(e) return( rep( 0, ncol(prior_ss50)) ) )
                                    } else out <- tryCatch( colMeans( exprMat[ match( targets , rownames(exprMat) ) , ], na.rm = T) , 
                                                            error=function(e) return(  rep( 0, ncol(prior_ss50))) )
                                  return( out )
                                    }) )
      rownames(TFA_mean) <- colnames(prior)
      return(TFA_mean)
    })
    
    a <- do.call(abind, c(TFA_average_50ss, list(along=3)))
    TFA_average_50ss <- apply(a, 1:2, median)
    
    # rownames(TFA_average_50ss) <- colnames(ATACseq_tgenesM_expr)
    # colnames(TFA_average_50ss) <- colnames(expMat_list[[ii]])
    
    return(TFA_average_50ss)
    gc()
  }
  
### calculate TFA with AUCell ####  
    TFAaucell_50ss_func <- function( exprMat,
                                  # geneSets= geneSets_ATACseq.priors ,
                                  TFs = TF_MeanMed ,
                                  ATACseq_prior = prior)
    {
      require(AUCell)
      require(GSEABase)
     

      # make gene sets  from TF targets
      genesets_list <- lapply( TFs , function(TF) {
        GeneSet( rownames( ATACseq_prior )[ ATACseq_prior[,TF] > 0 ], 
                 setName= TF ) 
        } )
      geneSets <-  GeneSetCollection( genesets_list )
    
       

        # 1. Build gene-expression rankings for each cell  
        cellRanks <- AUCell_buildRankings(  exprMat , nCores=1, plotStats=F)
        
        # 2. Calculate enrichment for the gene signatures (AUC)
        cells_AUC <- AUCell_calcAUC( geneSets , cellRanks )
        cells_AUC <- getAUC( cells_AUC)
        
        
        return(cells_AUC)
        gc()
      }


  
    


#### compute smothed mRNA signal and sort by PDS ####
TFA_MAW_func <- function( exprMatrices  , 
                          PDS_score  ,
                          wt_indeces  ,
                          WsizeF=20, WstepF=40 ,
                          expName , 
                          gnames )
  {
  ### function applies moving average window to smooth the signal
  ### and return 3 smooth objects: 1) nonscaled, 2) diff with wtype 3) scaled 
  # WsizeF - window size factor ( will define number of chunks  data will be devided )
  # WstepF - window step size factor 
  # exprMatrices - data of non smooth TFA, rows are genes, columns are cells
  # PDS_score - df with PDS score stored in a 'score' column 
  # and which rownames are same as colnames of the exprMatrices
  # gnames - TF that we want to investigate
  # expNames - experiment names (one per exprMatrix) for metadata
  
  # select TFs of interest
  exprMatrices <- exprMatrices[match( gnames , rownames(exprMatrices)),]
  
  ### non scaled moving average window
  {
    # add score column to the TFA matrix
    datTOplot <-  cbind.data.frame( t( exprMatrices )  , 
                                    score= PDS_score$score[ match( colnames( exprMatrices ) , 
                                                                   rownames( PDS_score) )] ) 
    # order cells by PDS
    datTOplot <- datTOplot[ order(datTOplot$score),]
    print(dim(datTOplot))
    
    # apply MAW to the df
    datTOplot <- as.data.frame ( (apply( datTOplot, 2 , function(x) slideFunct(x, 
                                                                               nrow(datTOplot)/WsizeF, 
                                                                               nrow(datTOplot)/WstepF ) ) ) )
    
  }
  
  ### to scale use the difference between the TFA lvl in a cell and 
  # the mean TFA lvl across all wtype cells of the same experiment
  {
    # calculate TFA diff for each cell
    exprMatrices_lfc <- apply(exprMatrices, 1, function( datt) {
      # get the mean expr of control cells
      meanWT <- mean( datt[ wt_indeces ] , na.rm = T ) 
      # print(meanWT)
      # XX <- log(datt/meanWT) 
      XX <- (datt-meanWT)
      XX[is.infinite(XX)] <- 0
      return(XX)}
    )
    
    # divide  by the average absolute difference of a given study
    # to make studies comparable
    exprMatrices_lfc <- exprMatrices_lfc/mean(abs(exprMatrices_lfc), na.rm = T)
    
    datTOplot_lfc <- cbind.data.frame( exprMatrices_lfc , 
                                       score=PDS_score$score[ 
                                         match( rownames( exprMatrices_lfc) , 
                                                rownames( PDS_score ) )] 
    )
    # sort by PDS
    datTOplot_lfc <- datTOplot_lfc[ order(datTOplot_lfc$score),]
    
    # use rolling average window to smooth the signal
    datTOplot_lfc <-  as.data.frame( apply( datTOplot_lfc, 2 ,
                                            function(x) slideFunct(x,
                                                                   nrow(datTOplot_lfc)/WsizeF,
                                                                   nrow(datTOplot_lfc)/WstepF ) 
    ) )
    
  }
  
  
  ### scale by rows
  {
    # scale the matrix
    datTOplot_scale <- scale( t( exprMatrices ) )
    datTOplot_scale <- cbind.data.frame(  datTOplot_scale  , 
                                          score=PDS_score$score[ 
                                            match( colnames( exprMatrices ) , 
                                                   rownames( PDS_score ) )] )
    
    datTOplot_scale <- datTOplot_scale[ order(datTOplot_scale$score),]
    
    datTOplot_scale <- as.data.frame ( apply( datTOplot_scale, 2 , 
                                              function(x) slideFunct(x, 
                                                                     nrow(datTOplot_scale)/WsizeF, 
                                                                     nrow(datTOplot_scale)/WstepF) 
    ) )
  }
  
  # assign same names for all columns
  colnames(datTOplot)[-ncol(datTOplot)] <- colnames(datTOplot_lfc)[-ncol(datTOplot_lfc)]  <- colnames(datTOplot_scale)[-ncol(datTOplot_scale)]  <-  gnames
  # add experiment metadata
  datTOplot$experiment <- datTOplot_lfc$experiment <- datTOplot_scale$experiment <- expName
  
  return( list( datTOplot , datTOplot_lfc, datTOplot_scale) )
}
