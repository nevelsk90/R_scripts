
# Wrapup that spits result table
RobertGO_wrapup <- function(geneset,IDconv=tx2geneGO,universe=universe)
{
  geneset <- unique(IDconv[IDconv[,1]%in%geneset,2])
  geneset <- as.character(geneset)
  geneset <- geneset[!is.na(geneset)]
  
  # If Universe is not global
  if (length(intersect(IDconv[,2],universe))==0) {universe <- as.character(universe) 
  universe <- unique(IDconv[IDconv[,1]%in%universe,2])
  universe <- geneset[!is.na(universe)]}
  
  RobertGO <- sf.clusterGoByGeneset(gomatrix,geneset,universe,min.genes=3)
  #RobertGO_q0.1prim=subset(RobertGO$results, Enrichment>0 & Is.primary==TRUE & Primary.Fisher.adj<0.1)[,1:13]
  # RobertGO_p0.05=subset(RobertGO$results,Fisher<0.05)[,1:13]
  return(RobertGO$results)
}

# string wrapper, to fit a list of Go terms on a panels
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



# FUNCTION
# shortens strings to n if they are longer than n
rs.truncateString = function( string, n=20, leeway=5, trail="...") {
  if (nchar(string) >= n ) {
    return( paste0( substr(string, 1, (n-leeway) ), trail ) )
  } else {
    return(string)
  }
}


# FUNCTION
# Based on BioConductor annotation pack, first construct a sparse matrix with gene members
# egGO2ALLEGS = annotation pack mapping of GO term to all associated gene identifiers
# Returns: Sparse binary matrix indicating membership of genes (columns) in GO terms (rows)
library("org.Mm.eg.db")
sf.createGoMatrix = function( egGO2ALLEGS = org.Mm.egGO2ALLEGS, min.genes=3, max.genes=1000, 
                              excludeIEA=F) {
  require(Matrix)
  
  golist = as.list(egGO2ALLEGS)
  golist = lapply( golist, function(x) {
    if (excludeIEA) { 
        this.iea = which( names(x) %in% "IEA" )
        if (length(this.iea) > 0 ) { x = x[-this.iea] }
      }
    return( unique(x) )
  })
  
  gn = sapply( golist,function(x) length( x ) )
  golist = golist[gn >= min.genes & gn <= max.genes]
  
  ugenes = unique( unlist(golist))
  
  gomat = Matrix(0, ncol=length(ugenes), nrow=length(golist), dimnames = list(names(golist), ugenes) )
  
  for(i in 1:nrow(gomat)) { gomat[i,] = colnames(gomat) %in% golist[[i]] }
  
  return(gomat)
}




# FUNCTION
# Then, subset the matrix based on the genes of interest
# Caveat: Specify minimum size of terms
sf.clusterGoByGeneset = function( gomatrix, geneset, universe, geneAlias=NULL, 
                                  min.genes=3, cut_max = 5, dynamic_cut=F,
                                  ontologies=c("BP","MF","CC"), two_sided=F,
                                  representative_by_silhouette = F) {
  require(AnnotationDbi)
  require(cluster)
  
  # for now, only support explicit universe
  universe = intersect(colnames(gomatrix), unique(universe) )
  gomatrix = gomatrix[,universe]
  
  # subset by ontologies
  gomatrix = gomatrix[ Ontology(rownames(gomatrix)) %in% ontologies, ]
  
  
  # make sure the geneset is congruent with the GO matrix (also, no duplicates!)
  geneset = unique(geneset)
  gidx = which(geneset %in% colnames(gomatrix))
  geneset = geneset[gidx]
  if(!is.null(geneAlias)) { geneAlias = geneAlias[gidx] }
  
  # subset according to min-genes (minimum refers to "significant" genes)
  n_in_term = apply(gomatrix, 1, sum)
  n_selected = apply(gomatrix[,geneset,drop=F], 1, sum)
  
  idx = which(n_selected >= min.genes)
  
  if (length(idx)==0) { warning("No GO terms with the specified min.genes were found for the subset of interest!"); return(NA) }
  
  gomatrix = gomatrix[idx,,drop=F]
  n_in_term = n_in_term[idx]
  n_selected = n_selected[idx]
  
  
  
  # Fisher tests for significant enrichment --> data frame with GO = rownames
  flist = list()
  for (k in 1:nrow(gomatrix)) {
    cmat = rbind( c(sum( gomatrix[k,] == 0 & !(colnames(gomatrix) %in% geneset)), sum(gomatrix[k,] == 1 & !(colnames(gomatrix) %in% geneset )) ),
                  c(sum( gomatrix[k,] == 0 & (colnames(gomatrix) %in% geneset) ), sum( gomatrix[k,] == 1 & (colnames(gomatrix) %in% geneset )) ) )
    
    expected = as.array(margin.table(cmat,1)) %*% t(as.array(margin.table(cmat,2))) / margin.table(cmat)
    enrichment = log2( cmat[2,2] / expected[2,2] )
    pval = fisher.test(cmat, alternative = ifelse(two_sided, "two.sided", "greater") )$p.value
    
    flist[k] = list(c("Significant"=cmat[2,1], 
                      "Expected"=round(expected[2,2], digits = 2) , 
                      "enrichment"=round(enrichment, digits = 2), 
                      "p.value"=pval))
  }
  flist = t( as.data.frame(flist) )
  rownames(flist) = rownames(gomatrix)
  
  
  # now subset the gomatrix to just the geneset of interest
  gomatrix = gomatrix[,geneset,drop=F]
  
  
  # CLUSTER AND SELECT PRIMARY TERMS (only if >1 go term selected)
  if (nrow(gomatrix) > 1) {
    # calculate manhattan distance and cluster
    godist = dist(gomatrix, method="manhattan")
    goclust = hclust(godist, method = "complete")
    
    
    # cut tree so that the edit distance between clustered GO terms is at most max.distance
    if (dynamic_cut) {
      ci = 1:cut_max
      maxindex = -Inf
      for (c in ci) {
        sindex = mean( as.matrix( silhouette( cutree(goclust, h = c + 0.5), godist ) )[,3] )
        if (sindex > maxindex) { maxindex = sindex; cut_max = c }
      }
    }
    
    #return( list(goclust, godist) )
    gopartition = cutree(goclust, h = cut_max + 0.5)
    
    gosil = as.matrix( silhouette( cutree(goclust, h = cut_max + 0.5), godist ) )
    
    if( length(gosil)>1 ) {
      rownames(gosil) = goclust$labels
      sindex = mean( as.matrix( silhouette( cutree(goclust, h = cut_max + 0.5), godist ) )[,3] )
    } else {
      sindex = NA
    }
    
    
    
    # Select as representative the smallest GO term in the cluster (by total annotations)
    if (!representative_by_silhouette) {
      reps = c()
      for ( c in 1:max(gopartition) ) {
        this = n_in_term[ names(gopartition)[which(gopartition == c)] ]
        
        if ( sum(this == min(this)) > 1 ) { 
          this = this[this == min(this)]
          this.rep = names(this)[which.min(nchar(Term(names(this))))]    # if several terms apply, choose the one that's... shortest
        } 
        else {
          this.rep = names(this)[which.min(this)]
        }
        reps = c(reps, this.rep)
      }
    } # else, select as representative the term with the highest silhouette index of its cluster
    else { reps = c()
    for ( c in 1:max(gopartition) ) {
      this = gosil[which(gosil[,"cluster"] == c),,drop=F]
      if ( sum(this[,"sil_width"] == max(this[,"sil_width"])) > 1 ) {
        this = this[this[,"sil_width"] == max(this[,"sil_width"]),]
        this.rep = rownames(this)[which.min(nchar(Term(rownames(this))))]    # if several terms apply, choose the one that's... shortest
      } 
      else {
        this.rep = rownames(this)[which.min(this[,"sil_width"])]
      }
      reps = c(reps, this.rep)
    }
    }
  } else {
    gopartition = setNames( 1, rownames(gomatrix) )
    reps = setNames( rownames(gomatrix), 1 )
  }
  
  
  # Summary data frame
  gm = data.frame( "Primary"=Term(reps[gopartition]),
                   "Is.primary" = names(n_in_term) %in% reps,
                   "Terms.in.cluster" = 0,
                   "GO.ID"=names(n_in_term),
                   "Term"=Term(names(n_in_term)), 
                   "Ontology"=Ontology(names(n_in_term)),
                   "Annotated"=n_in_term, 
                   "Significant"=n_selected,
                   "Expected"=flist[,"Expected"],
                   "Enrichment"=2^flist[,"enrichment"],
                   "log2Enrichment"=flist[,"enrichment"],
                   "Fisher"=flist[,"p.value"],
                   "Primary.Fisher.adj"=NA,
                   "Genes"="",
                   "Genes.in.secondary"="",
                   row.names = names(n_in_term), 
                   stringsAsFactors = F)
  gm$Primary.Fisher.adj[gm$Is.primary] = p.adjust( gm$Fisher[gm$Is.primary], method="BH" )
  
  # Add the names of significant genes in each term & number in cluster
  for (i in 1:nrow(gm)) {
    if( !is.null(geneAlias)) { gm$Genes[i] = list(geneAlias[ gomatrix[i,] == 1 ]) } else {
      gm$Genes[i] = list(geneset[ gomatrix[i,] == 1 ])
    }
    gm$Terms.in.cluster[i] = sum(gm$Primary == gm$Primary[i])
  }
  
  
  # For primary terms only, list genes that are included ONLY in secondary terms
  for (i in 1:nrow(gm)) {
    if( gm$Term[i] %in% gm$Primary ) {
      excl = setdiff( unique( unlist( gm[which(gm$Primary == gm$Term[i]),"Genes"] ) ),
                      gm[i,"Genes"][[1]] )
      if (length(excl) > 0) { gm$Genes.in.secondary[i] = list( excl ) }
    }
  }
  
  
  # For ease of plotting, translate the labels 
  if (nrow(gomatrix) > 1) { goclust$labels = paste0( strtrim( Term(goclust$labels), 30) ) } else {goclust = godist = sindex = NA}
  
  return( list("results"=gm, "background"=universe, "geneset"=geneset, "gomatrix"=gomatrix, "dists"=godist,
               "clustering"=goclust, "partition"=gopartition, "cut_at"=cut_max + 0.5, "silhouette_index"=sindex ) )
}


# FUNCTION
# Adds directional t-tests for significant genes per term
# geneID2entrez=biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id"), mart = mart)

sf.addDirectionality = function( resultframe, limmaframe, id="entrezgene_id", logfcid = "logFC" , tx2gene = geneID2entrez)
  {
  resultframe$t_up = 1
  resultframe$t_down = 1
  resultframe$meanLFC = NA
  
  x1=tx2gene[ match( rownames(limmaframe) , tx2gene$external_gene_name ) ,]
  # x1=x1[!duplicated(x1$ensembl_gene_id),]
  limmaframe=cbind(as.data.frame(limmaframe),x1)
  
  
  for (j in 1:nrow(resultframe)) {
    tmp = resultframe[j,"Genes"][[1]]
    #x1=tx2gene[tx2gene$ensembl_gene_id%in%rownames(limmaframe) | tx2gene$external_gene_name%in%rownames(limmaframe) ,c(2,5)]
# # older version (before 2020)
#     tmplfc = limmaframe[ which(limmaframe$entrezgene %in% tmp), 2]
    tmplfc = limmaframe[ which(limmaframe[,id] %in% tmp), logfcid]
    resultframe$t_up[j] = tryCatch( t.test(tmplfc, alternative = "greater")$p.value , error=function(e) 1 )
    resultframe$t_down[j] =  tryCatch( t.test(tmplfc, alternative = "less")$p.value ,error=function(e) 1 ) 
    resultframe$meanLFC[j] =  tryCatch( mean(tmplfc, na.rm=T) ,error=function(e) 0 )
  }
  
  
  prime = which( resultframe$Is.primary )
  resultframe$t_up.adj = NA
  resultframe$t_down.adj = NA
  resultframe$t_up.adj[prime] = p.adjust(resultframe$t_up[prime], method = "BH")
  resultframe$t_down.adj[prime] = p.adjust(resultframe$t_down[prime], method = "BH")
  
  return(resultframe)
}



# FUNCTION
# Wrapper for sf.clusterGoByGeneset, to be applied to a list of limma tables -- must have the appropriate mapping, e.g. entrez!
sf.clusterGoLimmaBatch = function( gomatrix, limmalist, id="entrez", id_pval="adj.P.Val", id_logfc = "logFC", 
                                   alpha_cutoff = 0.1, logFC_cutoff = 0, t_test=T, ...) {
  lout = lapply(limmalist, function(x) {
    vip = x[ which( x[,id_pval] <= alpha_cutoff & abs(x[,id_logfc]) >= logFC_cutoff ), id ]
    vip = vip[!is.na(vip)]
    bg = x[,id]
    bg = bg[!is.na(bg)]
    cat("|")
    
    if (length(vip) > 0) {
    return( sf.clusterGoByGeneset(gomatrix = gomatrix, 
                                  geneset = vip,
                                  universe = bg,
                                  ...) )
    } else {return(NA)}
  })
  
  if (t_test) {
    for (i in 1:length(lout)) {
      if(!is.na(lout[i])) {
        cat(".")
        lout[[i]]$results$t_up = 1
        lout[[i]]$results$t_down = 1
        lout[[i]]$results$meanLFC = NA
        
        for (j in 1:nrow(lout[[i]]$results)) {
          tmp = lout[[i]]$results[j,"Genes"][[1]]
          tmplfc = limmalist[[i]][ which(limmalist[[i]][,id] %in% tmp) ,"logFC" ]
          lout[[i]]$results$t_up[j] = t.test(tmplfc, alternative = "greater")$p.value
          lout[[i]]$results$t_down[j] = t.test(tmplfc, alternative = "less")$p.value
          lout[[i]]$results$meanLFC[j] = mean(tmplfc, na.rm=T)
        }
        
        lout[[i]]$results$t_up.adj = p.adjust(lout[[i]]$results$t_up, method = "BH")
        lout[[i]]$results$t_down.adj = p.adjust(lout[[i]]$results$t_down, method = "BH")
      }
    }
  }
  
  return(lout)
}






# FUNCTION
# allows plotting of multiple sfgo frames together, selecting only terms that are primary in at least one frame
sf.histogramWrapper = function( sflist, go.alpha=0.05 ) {
  common = unique( unlist( sapply( sflist, function(x) rownames(x)[x[,"Is.primary"] & x[,"Primary.Fisher.adj"] <= go.alpha ] ) ) )
  sflist = lapply(sflist, function(x) x[rownames(x) %in% common,] )
  return(sflist)
}




# Bonus: MDS map of terms
sf.plotEmbedding = function( godist, membership, perplexity=15 ) {
  
  cols = rainbow(length(unique(membership)))
  cols = cols[as.factor(membership)]
  
  require(Rtsne)
  md = Rtsne(godist, perplexity=perplexity)
  
  par(xpd=T, mar=c(5,5,5,18))
  plot(md$Y, col=cols, cex=1, pch=15)
  legend("topright", inset=c(-0.6,0), legend=unique(membership), fill=unique(cols), pch=15, title="Group")
  par(xpd=F)
  
  invisible(md)
}



# FUNCTION
margin_x=c(0,15,3,0)
# Plot a single result frame
sf.plotResults <- function( sfresult, alpha = 0.05, t.alpha = 0.05, max_bars = NULL, top_n = Inf,
                           only_positive_enrichment=T, order_by_enrichment=T, main="Enriched primary GO terms", 
                           legend.title = "DE genes", xlim=NULL, lab_max_chars=50, ... ) 
               {
  par(mar=par("mar") + margin_x, xpd=F)
  
  # Thresholding and sorting
  sfresult = sfresult[which(sfresult$Primary.Fisher.adj <= alpha),,drop=F ]
  if(only_positive_enrichment) { sfresult = sfresult[which(sfresult$log2Enrichment > 0),,drop=F ] }
  if(order_by_enrichment) { sfresult = sfresult[order(sfresult$log2Enrichment),,drop=F ] } else {
    sfresult = sfresult[order(sfresult$Fisher, decreasing = T),,drop=F ]
  }
  
  # top n if required
  if(top_n < nrow(sfresult)) {
    sfresult = sfresult[(nrow(sfresult)-top_n):nrow(sfresult),]
  }
  
  cols = rep( "gray", nrow(sfresult) )
  # If directionality has been determined via t-test, color bars accordingly
  if (!is.null(sfresult$t_up)) {
    cols[ sfresult$t_up.adj <= t.alpha ] = "coral"
    cols[ sfresult$t_down.adj <= t.alpha ] = "deepskyblue2"
  }
  
  # Add buffer bars to ensure consistency
  if ( !is.null(max_bars) ) { baradjust = max_bars - nrow(sfresult) } else { baradjust = 0 }
  
  # Plot the actual thing
  if(is.null(xlim)) { xlim = c(0, ceiling( max(sfresult$log2Enrichment) ) ) }
  bp = barplot( c( rep(NA, baradjust), sfresult$log2Enrichment ),  
                names.arg = c( rep(NA, baradjust), sapply( sfresult$Term, rs.truncateString, n=lab_max_chars ) ), 
                horiz = T, las=1, 
                col = c( rep(NA, baradjust), "deepskyblue"), 
                lwd=2,
                xlim = xlim, 
                axes = F, ... )
  axis(side = 3)
  mtext("log2(enrichment)", line=2.5)
  title(main, line = 5, ...)
  
  # And add a legend
  if (!is.null(sfresult$t_up)) {
    legend(x = max(xlim)-3, y=bp[baradjust+1], fill=c("coral","deepskyblue2","gray23"),
           legend=c("up","down","indet."), xpd = T, bty = "n", y.intersp = 0.9, title=legend.title)
  }
  
  # add number of genes
  par(xpd=T)
  text(x = c( rep(NA, baradjust), sfresult$log2Enrichment ), y = bp, 
       c(rep(NA, baradjust), paste0("*", sfresult$Significant, "/", sfresult$Annotated) ), pos = 4)
  
  par(mar=par("mar") - margin_x, xpd=F)
  
  invisible(sfresult)
}


# FUNCTION
# plots mean logFC for terms that are significantly biased towards up-/down-regulation
sf.plotDirectionality = function( results, alpha=0.05, max_bars=NULL, Nterms=NULL, primeOnly=F, marginsPar=c(0,16,5,4),... ) 
  {
  if ( !is.null(Nterms) ) {
    results <- results[ which( results$Is.primary ), ,drop=F] 
    results <- head(results[order(results$Primary.Fisher.adj),], Nterms)
    } else results = results[ which( (results$t_up.adj <= alpha | results$t_down.adj <= alpha) & results$Is.primary ), ,drop=F]
  print(nrow(results))
  

  if (nrow(results) >= 1) {
    
    results = results[order(results$meanLFC),]
    cols = rep("deepskyblue2",nrow(results))
    cols[ results$meanLFC >= 0 ] = "coral"
    
    # Add buffer bars to ensure consistency
    if ( !is.null(max_bars) ) { baradjust = max_bars - nrow(results) } else { baradjust = 0 }
    
    
    par(mar=par("mar") + marginsPar, xpd=F)
    barplot(c( rep(NA, baradjust), results$meanLFC),
            names.arg = c( rep(NA, baradjust), results$Primary), 
            col = c( rep(NA,baradjust), cols),
            horiz = T, las=1, axes = F, ...)
    axis(side = 3)
    mtext("mean logFC of differentially regulated", line=2.5)
    abline(v=0)
    par(mar=par("mar") - marginsPar , xpd=F)
    title(..., line=2)
  } else {
    warning("No significantly regulated terms found!")
  }
}





# FUNCTION
# Beeswarm plots for GO terms
sf.plotSingleTermStripchart = function( gomatrix=NULL, goid=NULL, expressionframe, id_logFC = "logFC", id_pval = "adj.P.Val", id="entrez", alpha=0.1,
                                        lim=NULL, draw.boxplot=F, ylab="logFC", highlight=NULL, highlight.col="black", highlight.legend = "Highlighted",
                                        max.title.chars=42, title="") 
  {
  require("beeswarm")
  
  # Select genes belonging to term, and significant ones
  if (!is.null(gomatrix)) {
    vip = colnames(gomatrix)[ which(gomatrix[goid,] == 1) ] %>% intersect(expressionframe[,id])
    vip_significant = intersect(expressionframe[ which( expressionframe[,id_pval] <= alpha) ,id ], vip)
  } else {
    vip = expressionframe[,id]
    vip_significant = intersect(expressionframe[ which( expressionframe[,id_pval] <= alpha) ,id ], vip)
  }
  
  # select logFC's, sorted into sig/non-sig
    lfc_nonsig = setNames( expressionframe[ which( expressionframe[,id] %in% setdiff( vip, vip_significant) ), id_logFC ],
                           expressionframe[ which( expressionframe[,id] %in% setdiff( vip, vip_significant) ), id ] )
    lfc_sig = setNames( expressionframe[ which( expressionframe[,id] %in% vip_significant ), id_logFC ],
                        expressionframe[ which( expressionframe[,id] %in% vip_significant ), id ] )
    lfc_all = c( lfc_nonsig, lfc_sig )
  
  # set basic colors for up/down/neither
    lfc_scols = rep("coral",length(lfc_sig))
    lfc_scols[lfc_sig < 0] = "deepskyblue2"
    lfc_cols = c( rep("grey", length(lfc_nonsig)), lfc_scols)
    
  # adjust plot range limit
    if (is.null(lim)) { lim = c(-1,1) * max(abs(c(lfc_nonsig, lfc_sig)), na.rm = T) }
  
  # For points outside the limit, select special symbols to indicate that
    lfc_ub = which( lfc_all > max(lim, na.rm = T) )
    lfc_lb = which( lfc_all < min(lim, na.rm = T) )
    lfc_pchs = setNames( rep(20, length(lfc_all)) , names(lfc_all) )
    lfc_pchs[lfc_ub] = 2
    lfc_pchs[lfc_lb] = 6
  
  
  # Adjust symbols for highlighted genes, if any
    if( !is.null(highlight) ) {
      lfc_pchs[ lfc_pchs == 20 & names(lfc_pchs) %in% highlight ] = 21
      lfc_pchs[ lfc_pchs == 2 & names(lfc_pchs) %in% highlight ] = 24
      lfc_pchs[ lfc_pchs == 6 & names(lfc_pchs) %in% highlight ] = 25
    }
  
  
  # draw an additional boxplot if desired
    if(draw.boxplot) { boxplot( lfc_all, ylim = lim, border="black", outline=FALSE ) }
  
  # set values outside the limit to the closest limit
    lfc_all[ lfc_all > max(lim, na.rm = T) ] = max(lim, na.rm = T)
    lfc_all[ lfc_all < min(lim, na.rm = T) ] = min(lim, na.rm = T)
    
    lfc_all[is.na(lfc_all)] = 0
    
    
    
  # make the beeswarm plot
    beeswarm( lfc_all, 
              pwcol = lfc_cols, bg=highlight.col, 
              pwpch = lfc_pchs,
              add=draw.boxplot,
              ylim=lim, axes=F, corral = "gutter")
    abline(h=0, lty=2, col="grey")
    box(col="grey")
    axis(side = 4, las=1,  mgp=c(3, 0.5, 0), tck=-0.03, lwd=2)
    
  
  # title and axis label
    tmp = Term(goid)
    if( nchar(tmp) > max.title.chars ) {
      tmp = paste0( substr(tmp, 1, max.title.chars - 3), "... ")
    }
    mtext(side = 2, paste0(tmp," (", goid,")"), line = 0.2, col="dodgerblue4")
    mtext(side = 4, ylab, line = 1.3)
  
  # text: Significant / total
    mtext(side = 1, text = paste0("*", length(lfc_sig),"/", length(lfc_all) ) )
    
  # legend for highlight
    if (!is.null(highlight)) {
      
      legend("bottomleft", pch=c(19,19,19,19), col=c("deepskyblue2","coral","black","grey"), 
             legend=c("up-regulated","down-regulated","Wol-dependent","not sig."), bty="n")
    }
    
  # Extra title
    title(title)
    
    
    invisible( expressionframe[ expressionframe[,id] %in% names(lfc_all), ])
}

# FUNCTION
# make beeswarm plots for every significantly differential primary term in goresults
sf.plotBatchStripchart = function( gomatrix, goresults, expressionframe, id_logFC = "logFC", id_pval = "adj.P.Val", id="entrez", alpha=0.1, lim=NULL,
                                   directionality_alpha = 0.05, highlight = NULL, ylab="logFC", highlight.legend = "Highlighted", title="" ) {
  select_terms = rownames(goresults)[ ( goresults$t_up.adj <= directionality_alpha | goresults$t_down.adj <= directionality_alpha ) & goresults$Is.primary ]
  select_terms = select_terms[!is.na(select_terms)]
  
  for (t in select_terms) {
    print(t)
    sf.plotSingleTermStripchart(gomatrix = gomatrix, goid = t, lim=lim,
                                expressionframe = expressionframe,
                                id_logFC = id_logFC,
                                id_pval = id_pval, 
                                id = id, 
                                alpha = alpha, 
                                highlight = highlight,
                                highlight.legend = highlight.legend,
                                ylab = ylab,
                                title=title)
  }
}


# FUNCTION
# helper function
add.alpha = function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


# Bonus: Kappa statistic as in DAVID
# Conclusion: Gets biased by the fact that many genes are not included in ANY terms of any given cluster... seems like an oversight, at DAVID, shocking!
rs.kappa = function( v1, v2 ) {
  p0 = sum(v1 == v2) / length(v1)
  ma = ( ( sum( v1 == 1 & v1 == v2 ) + sum( v1 == 0 & v1 != v2 ) ) * ( sum( v1 == 1 & v1 == v2 ) + sum( v2 == 0 & v1 != v2 ) ) ) / length(v1)
  mb = ( ( sum( v2 == 1 & v2 == v1 ) + sum( v2 == 0 & v2 != v1 ) ) * ( sum( v2 == 1 & v2 == v1 ) + sum( v1 == 0 & v2 != v1 ) ) ) / length(v2)
  pe = (ma + mb) / length(v1)
  kappa = (p0 - pe) / (1-pe)
  return(kappa)
}

# calculate kappa distances between rows, just like dist
rs.kappa.dist = function(m) {
  md = matrix(0, nrow = nrow(m), ncol = nrow(m))
  for (i in 1:nrow(md)) {
    for (j in i:ncol(md)) {
      md[j,i] = rs.kappa(m[i,], m[j,])
    }
  }
  return(as.dist(md))
}












