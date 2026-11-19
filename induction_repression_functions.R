get_induction_repression <- function(second_derivate=data.frame(),
                                     gene_names_cutoff=character(),
                                     slopes_frame_TC=data.frame(),
                                     MinMaxFrame=data.frame(),
                                     MinMaxType="log2FC"){
  
  ## function will calculate induction and repression time points for time series based on
  ## slopes and changes of slopes over the time series for a set of genes
  
  ## The function does miss genes that are strictly monotonic in their slopes if they are changing
  ## more or less constant over the time. You can define genes that behave strictly monotonic but
  ## and change over time but are not assigned induction or repression time points
  ##and choose 0 as induction or repression time point for those
  
  ## second_derivate should be a dataframe that has the changes in slopes of neighbouring time points
  ## as columns and genes in rows
  
  ## slopes_frame_TC should be a dataframe that has the slopes at given time points
  ## as columns and genes in rows
  ## colnames should have structure: slope_timepoint1_timepoint2
  
  ## number of colums in second_derivate should be the same as in slopes_frame_TC
  ## an additional column will be added to the beginning of slopes_frame_TC (slope_0_0)
  ## with the assumption that the gene is coming out of a steady state
  
  ## gene_names_cutoff is a character vector with names of genes that should be analyzed as the
  ## function will identify mainly noise if genes do not change over the time course and induction
  ## and repression points are not usefull
  
  ## MinMaxFrame is the data frame that was used to calculate slopes in the first place
  ## might be TPMs or log2FCs at the moment (important to define a dynamic threshold per gene)
  
  
  ## in this step the 2nd derivative and slope frames are cut to the genes that are defined by
  ## gene_names_cutoff and a slope_0_0 column is added to the beginning of the slopes frame
  holder <- trim_2nd_derivative_and_slopes(second_derivate,
                                           gene_names_cutoff,
                                           slopes_frame_TC)
  
  short_2nd_deriv <- holder[[1]]
  short_slopes_frame_TC <- holder[[2]]
  MinMaxFrame <- MinMaxFrame[which(row.names(MinMaxFrame)%in%
                                     gene_names_cutoff),]
  
  ## which columns in the second derivate frame exceed a dynamic threshold and migh be considered
  ## as an induction or repression time point
  cols_2nd_deriv <- list()
  
  for(d in 1:dim(short_2nd_deriv)[1]){
    gene <- row.names(short_2nd_deriv)[d]
    
    if(MinMaxType=="TPM"){
      dynamic_thresh <- log2(((max(MinMaxFrame[gene,])+0.01)/(min(MinMaxFrame[gene,])+0.01)))/(dim(short_slopes_frame_TC)[2]-1)
    }
    else if(MinMaxType=="log2FC"){
      dynamic_thresh <- dist(c(max(MinMaxFrame[gene,]),min(MinMaxFrame[gene,])))/(dim(short_slopes_frame_TC)[2]-1)
    }
    
    holder <- short_2nd_deriv[gene,]
    holder <- holder[which(abs(holder) > dynamic_thresh)]
    names(holder) <- which(abs(short_2nd_deriv[gene,]) > dynamic_thresh)
    
    if(length(holder)>0){
      cols_2nd_deriv <- c(cols_2nd_deriv, list(holder))
    } else {
      cols_2nd_deriv <- c(cols_2nd_deriv, list(NA))
    }
    
  }
  names(cols_2nd_deriv) <- row.names(short_2nd_deriv)
  
  ## list of cols that might be induction or repression time points and slope and 2nd
  ## derivative  frames are used to identify which of those points are actually an
  ## induction or repression point
  induction_repression_all_points <- get_induction_repression_points(cols_2nd_deriv,
                                                                     short_2nd_deriv,
                                                                     short_slopes_frame_TC)
  
  return(induction_repression_all_points)
}

trim_2nd_derivative_and_slopes <- function(second_derivate=data.frame(),
                                           gene_names_cutoff=character(),
                                           slopes_frame_TC=data.frame()){
  
  ## only selecting genes that were predefined (genes that change over time)
  ## in the 2nd derivative and slope frames
  ## also adding a slope_0_0 column to the slope frame
  short_2nd_deriv <- second_derivate[which(row.names(second_derivate)%in%
                                             gene_names_cutoff),]
  
  short_slopes_frame_TC <- slopes_frame_TC[which(row.names(slopes_frame_TC)%in%
                                                   gene_names_cutoff),]
  short_slopes_frame_TC <- cbind(matrix(0, nrow = dim(short_slopes_frame_TC)[1],
                                        ncol = 1), short_slopes_frame_TC)
  colnames(short_slopes_frame_TC)[1] <- "slope_0_0"
  
  return(list(short_2nd_deriv,short_slopes_frame_TC))
}

get_induction_repression_points <- function(cols_2nd_deriv=data.frame(),
                                            short_2nd_deriv=data.frame(),
                                            short_slopes_frame_TC=data.frame()){
  
  ## function that uses information of slopes and 2nd derivative from the time course as well
  ## as information on which columns surpass a certain change in slopes (2nd derivative).
  ## the output will be a dataframe that is ready to be plotted with ggplot. Each row is
  ## one induction or repression time point which has information on the type (induction or repression),
  ## the Gene (row names of dataframes) and the time point (based on col names of slopes frame)
  
  induction_repression_all_points <- list()
  for(d in 1:length(cols_2nd_deriv)){
    if(!is.na(cols_2nd_deriv[[d]])[1]){
      holder <- c()
      for(c in 1:length(cols_2nd_deriv[[d]])){
        if(short_slopes_frame_TC[d,(as.numeric(names(cols_2nd_deriv[[d]]))[c]+1)]<0){
          holder <- c(holder, -1)
        }
        else if(short_slopes_frame_TC[d,(as.numeric(names(cols_2nd_deriv[[d]])[c]))+1]>0){
          holder <- c(holder, 1)
        } else {
          holder <- c(holder, 0)
        }
      }
      names(holder) <- sapply(strsplit(colnames(short_slopes_frame_TC)[as.numeric(names(cols_2nd_deriv[[d]]))+1],
                                       split = "_"),function(x) x[2])
      
      ## we only want to define the first time point to be an induction or repression time point
      ## if it is a series of time points that make the cut off. We also don't want an induction
      ## point followed by an induction point (same for repression)
      omit <- c()
      if(length(holder)>1){
        for(p in 2:length(holder)){
          if(holder[p]==holder[p-1]){
            omit <- c(omit,p)
          }
        }
        holder <- holder[which(!1:length(holder)%in%omit)]
      }
      
      induction_repression_all_points <- c(induction_repression_all_points,
                                           list(holder))
    } else {
      induction_repression_all_points <- c(induction_repression_all_points,
                                           list(NA))
    }
  }
  names(induction_repression_all_points) <- row.names(short_2nd_deriv)
  
  return(induction_repression_all_points)
}