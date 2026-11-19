### function to convert long df in a wide df for GO chord plot
chord_datX <- function(data, genes, process){
  # require(GOplot)
  
  id <- term <- logFC <- BPprocess <- NULL
  
  # colnames(data) <- tolower(colnames(data))
  if (missing(genes)){
    if (is.null(data$logFC)){
      genes <- as.character(unique(data$genes))
    }else{
      genes <- subset(data, !duplicated(genes), c(genes, logFC))
    }
  }else{
    if(is.vector(genes)){
      genes <- as.character(genes) 
    }else{
      if(class(genes[, 2]) != 'numeric') genes[, 2] <- as.numeric(levels(genes[, 2]))[genes[, 2]]
      genes[, 1] <- as.character(genes[, 1])
      colnames(genes) <- c('genes', 'logFC')
    }
  }
  if (missing(process)){
    process <- as.character(unique(data$term))
  }else{
    if(class(process) != 'character') process <- as.character(process)
  }
  if (strsplit(process[1],':')[[1]][1] == 'GO'){
    subData <- subset(data, id%in%process)
    colnames(subData)[which(colnames(subData) == 'id')] <- 'BPprocess'
  }else{
    subData <- subset(data, term%in%process)
    colnames(subData)[which(colnames(subData) == 'term')] <- 'BPprocess'
  }
  
  if(is.vector(genes)){
    M <- genes[genes%in%unique(subData$genes)]
    mat <- matrix(0, ncol = length(process), nrow = length(M))
    rownames(mat) <- M
    colnames(mat) <- process
    for (p in 1:length(process)){
      sub2 <- subset(subData, BPprocess == process[p])
      for (g in 1:length(M)) mat[g, p] <- ifelse(M[g]%in%sub2$genes, 1, 0)
    }
  }else{
    genes <- subset(genes, genes %in% unique(subData$genes))
    N <- length(process) + 1
    M <- genes[,1] 
    mat <- matrix(0, ncol = N, nrow = length(M))
    rownames(mat) <- M
    colnames(mat) <- c(process, 'logFC') 
    mat[,N] <- genes[,2]
    for (p in 1:(N-1)){
      sub2 <- subset(subData, BPprocess == process[p])
      for (g in 1:length(M)) mat[g, p] <- ifelse(M[g]%in%sub2$genes, 1, 0)
    }
  }
  return(mat)
}

### MAIN circular chord-plot function
GOChordX <- function(data, title, space, gene.order, gene.size, gene.space, nlfc = 1, lfc.col, lfc.min, lfc.max, ribbon.col, border.size, process.label, limit)
{
  require(scales)
  #require(GOplot)
  y <- id <- xpro <- ypro <- xgen <- ygen <- lx <- ly <- ID <- logFC <- NULL
  Ncol <- dim(data)[2]
  
  if (missing(title)) title <- ''
  if (missing(space)) space = 0
  if (missing(gene.order)) gene.order <- 'none'
  if (missing(gene.size)) gene.size <- 3
  if (missing(gene.space)) gene.space <- 0.2
  if (missing(lfc.col)) lfc.col <- c('blue', 'white', 'red')
  if (missing(lfc.min)) lfc.min <- min(data[,"logFC"], na.rm = T)
  if (missing(lfc.max)) lfc.max <-  max(data[,"logFC"], na.rm = T)
  if (missing(border.size)) border.size <- 0.5
  if (missing (process.label)) process.label <- 11
  if (missing(limit)) limit <- c(0, 0)
  
  if (gene.order == 'logFC') data <- data[order(data[, Ncol], decreasing = T), ]
  if (gene.order == 'alphabetical') data <- data[order(rownames(data)), ]
  if (sum(!is.na(match(colnames(data), 'logFC'))) > 0){
    if (nlfc == 1){
      cdata <- check_chord(data[, 1:(Ncol - 1)], limit)
      lfc <- sapply(rownames(cdata), function(x) data[match(x,rownames(data)), Ncol])
    }else{
      cdata <- check_chord(data[, 1:(Ncol - nlfc)], limit)
      lfc <- sapply(rownames(cdata), function(x) data[, (Ncol - nlfc + 1)])
    }
  }else{
    cdata <- check_chord(data, limit)
    lfc <- 0
  }
  
  if (missing(ribbon.col)) colRib <- grDevices::rainbow(dim(cdata)[2]) else colRib <- ribbon.col
  nrib <- colSums(cdata)
  ngen <- rowSums(cdata)
  Ncol <- dim(cdata)[2]
  Nrow <- dim(cdata)[1]
  colRibb <- c()
  for (b in 1:length(nrib)) colRibb <- c(colRibb, rep(colRib[b], 202 * nrib[b]))
  r1 <- 1; r2 <- r1 + 0.1
  xmax <- c(); x <- 0
  for (r in 1:length(nrib)){
    perc <- nrib[r] / sum(nrib)
    xmax <- c(xmax, (pi * perc) - space)
    if (length(x) <= Ncol - 1) x <- c(x, x[r] + pi * perc)
  }
  
  xp <- c(); yp <- c()
  l <- 50
  
  for (s in 1:Ncol){
    xh <- seq(x[s], x[s] + xmax[s], length = l)
    xp <- c(xp, r1 * sin(x[s]), r1 * sin(xh), r1 * sin(x[s] + xmax[s]), r2 * sin(x[s] + xmax[s]), r2 * sin(rev(xh)), r2 * sin(x[s]))
    yp <- c(yp, r1 * cos(x[s]), r1 * cos(xh), r1 * cos(x[s] + xmax[s]), r2 * cos(x[s] + xmax[s]), r2 * cos(rev(xh)), r2 * cos(x[s]))
  }
  
  df_process <- data.frame(x = xp, y = yp, id = rep(c(1:Ncol), each = 4 + 2 * l))
  xp <- c(); yp <- c(); logs <- NULL
  x2 <- seq(0 - space, -pi - (-pi / Nrow) - space, length = Nrow)
  xmax2 <- rep(-pi / Nrow + space, length = Nrow)
  
  for (s in 1:Nrow){
    xh <- seq(x2[s], x2[s] + xmax2[s], length = l)
    if (nlfc <= 1){
      xp <- c(xp, (r1 + 0.05) * sin(x2[s]), (r1 + 0.05) * sin(xh), (r1 + 0.05) * sin(x2[s] + xmax2[s]), r2 * sin(x2[s] + xmax2[s]), r2 * sin(rev(xh)), r2 * sin(x2[s]))
      yp <- c(yp, (r1 + 0.05) * cos(x2[s]), (r1 + 0.05) * cos(xh), (r1 + 0.05) * cos(x2[s] + xmax2[s]), r2 * cos(x2[s] + xmax2[s]), r2 * cos(rev(xh)), r2 * cos(x2[s]))
    }else{
      tmp <- seq(r1, r2, length = nlfc + 1)
      for (t in 1:nlfc){
        logs <- c(logs, data[s, (dim(data)[2] + 1 - t)])
        xp <- c(xp, (tmp[t]) * sin(x2[s]), (tmp[t]) * sin(xh), (tmp[t]) * sin(x2[s] + xmax2[s]), tmp[t + 1] * sin(x2[s] + xmax2[s]), tmp[t + 1] * sin(rev(xh)), tmp[t + 1] * sin(x2[s]))
        yp <- c(yp, (tmp[t]) * cos(x2[s]), (tmp[t]) * cos(xh), (tmp[t]) * cos(x2[s] + xmax2[s]), tmp[t + 1] * cos(x2[s] + xmax2[s]), tmp[t + 1] * cos(rev(xh)), tmp[t + 1] * cos(x2[s]))
      }}}
  if(lfc[1] != 0){
    if (nlfc == 1){
      df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow), each = 4 + 2 * l), logFC = rep(lfc, each = 4 + 2 * l))
    }else{
      df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:(nlfc*Nrow)), each = 4 + 2 * l), logFC = rep(logs, each = 4 + 2 * l))  
    }
  }else{
    df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow), each = 4 + 2 * l))
  }
  aseq <- seq(0, 180, length = length(x2)); angle <- c()
  for (o in aseq) if((o + 270) <= 360) angle <- c(angle, o + 270) else angle <- c(angle, o - 90)
  df_texg <- data.frame(xgen = (r1 + gene.space) * sin(x2 + xmax2/2),ygen = (r1 + gene.space) * cos(x2 + xmax2 / 2),labels = rownames(cdata), angle = angle)
  df_texp <- data.frame(xpro = (r1 + 0.15) * sin(x + xmax / 2),ypro = (r1 + 0.15) * cos(x + xmax / 2), labels = colnames(cdata), stringsAsFactors = FALSE)
  cols <- rep(colRib, each = 4 + 2 * l)
  x.end <- c(); y.end <- c(); processID <- c()
  for (gs in 1:length(x2)){
    val <- seq(x2[gs], x2[gs] + xmax2[gs], length = ngen[gs] + 1)
    pros <- which((cdata[gs, ] != 0) == T)
    for (v in 1:(length(val) - 1)){
      x.end <- c(x.end, sin(val[v]), sin(val[v + 1]))
      y.end <- c(y.end, cos(val[v]), cos(val[v + 1]))
      processID <- c(processID, rep(pros[v], 2))
    }
  }
  df_bezier <- data.frame(x.end = x.end, y.end = y.end, processID = processID)
  df_bezier <- df_bezier[order(df_bezier$processID,-df_bezier$y.end),]
  x.start <- c(); y.start <- c()
  for (rs in 1:length(x)){
    val<-seq(x[rs], x[rs] + xmax[rs], length = nrib[rs] + 1)
    for (v in 1:(length(val) - 1)){
      x.start <- c(x.start, sin(val[v]), sin(val[v + 1]))
      y.start <- c(y.start, cos(val[v]), cos(val[v + 1]))
    }
  }	
  df_bezier$x.start <- x.start
  df_bezier$y.start <- y.start
  df_path <- bezier(df_bezier, colRib)
  # if(length(df_genes$logFC) != 0){
  #   tmp <- sapply(df_genes$logFC, function(x) ifelse(x > lfc.max, lfc.max, x))
  #   logFC <- sapply(tmp, function(x) ifelse(x < lfc.min, lfc.min, x))
  #   df_genes$logFC <- logFC
  # }
  
  g<- ggplot() +
    geom_polygon(data = df_process, aes(x, y, group=id), fill='gray70', inherit.aes = F,color='black') +
    geom_polygon(data = df_process, aes(x, y, group=id), fill=cols, inherit.aes = F,alpha=0.6,color='black') +	
    geom_point(aes(x = xpro, y = ypro, size = factor(labels, levels = labels), shape = NA), data = df_texp) +
    guides(size = guide_legend("GO Terms", ncol = 2 , byrow = T, override.aes = list(shape = 22, fill = unique(cols), size = 8))) +
    theme(legend.text = element_text(size = process.label)) +
    geom_text(aes(xgen, ygen, label = labels, angle = angle), data = df_texg, size = gene.size) +
    geom_polygon(aes(x = lx, y = ly, group = ID), data = df_path, fill = colRibb, color = 'black', size = border.size, inherit.aes = F) +		
    labs(title = title) +
    theme_blank
  
  if (nlfc >= 1){
    g + geom_polygon(data = df_genes, aes(x, y, group = id, fill = logFC), inherit.aes = F, color = 'black') +
      # scale_fill_gradient2('logFC', space = 'Lab', low = lfc.col[3], mid = lfc.col[2], high = lfc.col[1], guide = guide_colorbar(title.position = "top", title.hjust = 0.5),
      #                      breaks = c(min(df_genes$logFC), max(df_genes$logFC)), labels = c(round(min(df_genes$logFC)), round(max(df_genes$logFC)))) +
      scale_fill_gradientn(colours = c("blue","white","red"), values = rescale( c(lfc.min, 0 , lfc.max) ),
                           guide = "colorbar", name='log2FC' ,limits=c( lfc.min, lfc.max) ) +  # color for pcorr as corr
      theme(legend.position = 'bottom', legend.background = element_rect(fill = 'transparent'), legend.box = 'horizontal', legend.direction = 'horizontal')
  }else{
    g + geom_polygon(data = df_genes, aes(x, y, group = id), fill = 'gray50', inherit.aes = F, color = 'black')+
      theme(legend.position = 'bottom', legend.background = element_rect(fill = 'transparent'), legend.box = 'horizontal', legend.direction = 'horizontal')
  }
}

##############
# In general #
##############

# Theme blank
theme_blank <- theme(axis.line = element_blank(), axis.text.x = element_blank(),
                     axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
                     axis.title.y = element_blank(), panel.background = element_blank(), panel.border = element_blank(),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())

# Draw adjacent table for GOBubble and GOCircle
draw_table <- function(data, col){
  id <- term <- NULL
  colnames(data) <- tolower(colnames(data))
  if (missing(col)){
    tt1 <- ttheme_default()
  }else{
    text.col <- c(rep(col[1], sum(data$category == 'BP')), rep(col[2], sum(data$category == 'CC')), rep(col[3], sum(data$category == 'MF')))
    tt1 <- ttheme_minimal(
      core = list(bg_params = list(fill = text.col, col=NA, alpha= 1/3)), 
      colhead = list(fg_params = list(col = "black")))
  }
  table <- tableGrob(subset(data, select = c(id, term)), cols = c('ID', 'Description'), rows = NULL, theme = tt1)
  return(table)
}

###########
# GOChord #
###########

# Bezier function for drawing ribbons
bezier <- function(data, process.col){
  x <- c()
  y <- c()
  Id <- c()
  sequ <- seq(0, 1, by = 0.01)
  N <- dim(data)[1]
  sN <- seq(1, N, by = 2)
  if (process.col[1] == '') col_rain <- grDevices::rainbow(N) else col_rain <- process.col
  for (n in sN){
    xval <- c(); xval2 <- c(); yval <- c(); yval2 <- c()
    for (t in sequ){
      xva <- (1 - t) * (1 - t) * data$x.start[n] + t * t * data$x.end[n]
      xval <- c(xval, xva)
      xva2 <- (1 - t) * (1 - t) * data$x.start[n + 1] + t * t * data$x.end[n + 1]
      xval2 <- c(xval2, xva2)
      yva <- (1 - t) * (1 - t) * data$y.start[n] + t * t * data$y.end[n]  
      yval <- c(yval, yva)
      yva2 <- (1 - t) * (1 - t) * data$y.start[n + 1] + t * t * data$y.end[n + 1]
      yval2 <- c(yval2, yva2)			
    }
    x <- c(x, xval, rev(xval2))
    y <- c(y, yval, rev(yval2))
    Id <- c(Id, rep(n, 2 * length(sequ)))
  }
  df <- data.frame(lx = x, ly = y, ID = Id)
  return(df)
}

# Check function for GOChord argument 'limit'
check_chord <- function(mat, limit){
  
  if(all(colSums(mat) >= limit[2]) & all(rowSums(mat) >= limit[1])) return(mat)
  
  tmp <- mat[(rowSums(mat) >= limit[1]),]
  mat <- tmp[,(colSums(tmp) >= limit[2])]
  
  mat <- check_chord(mat, limit)
  return(mat)
}

##########
# GOVenn #
##########

# Calculate points to draw a circle
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

# Calculate overlap for three lists
get_overlap<-function(A,B,C){
  colnames(A)<-c('ID','logFC')
  colnames(B)<-c('ID','logFC')
  colnames(C)<-c('ID','logFC')
  UP<-NULL;DOWN<-NULL;Change<-NULL
  if (class(A$logFC)!='numeric'){
    A$logFC<-gsub(",", ".", gsub("\\.", "", A$logFC))
    A$Trend<-sapply(as.numeric(A$logFC), function(x) ifelse(x > 0,'UP','DOWN')) 
  }else{ A$Trend<-sapply(A$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
  if (class(B$logFC)!='numeric'){
    B$logFC<-gsub(",", ".", gsub("\\.", "", B$logFC))
    B$Trend<-sapply(as.numeric(B$logFC), function(x) ifelse(x > 0,'UP','DOWN')) 
  }else{ B$Trend<-sapply(B$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
  if (class(C$logFC)!='numeric'){
    C$logFC<-gsub(",", ".", gsub("\\.", "", C$logFC))
    C$Trend<-sapply(as.numeric(C$logFC), function(x) ifelse(x > 0,'UP','DOWN')) 
  }else{ C$Trend<-sapply(C$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
  if (sum(((A$ID%in%B$ID)==T)==T)==0){
    AB<-data.frame() 
  }else{
    AB<-A[(A$ID%in%B$ID)==T,which(colnames(A)%in%c('ID','logFC','Trend'))]
    BA<-B[(B$ID%in%A$ID)==T,which(colnames(B)%in%c('ID','logFC','Trend'))]
    AB<-merge(AB,BA,by="ID")
    rownames(AB)<-AB$ID
    AB<-AB[,-1]   
  }
  if (sum(((A$ID%in%C$ID)==T)==T)==0){
    AC<-data.frame() 
  }else{
    AC<-A[(A$ID%in%C$ID)==T,which(colnames(A)%in%c('ID','logFC','Trend'))]
    CA<-C[(C$ID%in%A$ID)==T,which(colnames(C)%in%c('ID','logFC','Trend'))]
    AC<-merge(AC,CA,by="ID")
    rownames(AC)<-AC$ID
    AC<-AC[,-1]
  }
  if (sum(((B$ID%in%C$ID)==T)==T)==0){
    BC<-data.frame() 
  }else{
    BC<-B[(B$ID%in%C$ID)==T,which(colnames(B)%in%c('ID','logFC','Trend'))]
    CB<-C[(C$ID%in%B$ID)==T,which(colnames(C)%in%c('ID','logFC','Trend'))]
    BC<-merge(BC,CB,by="ID")
    rownames(BC)<-BC$ID
    BC<-BC[,-1]
  }
  if (sum(((A$ID%in%B$ID)==T & (A$ID%in%C$ID)==T))==0){
    ABC<-data.frame() 
  }else{
    ABC<-A[((A$ID%in%B$ID)==T & (A$ID%in%C$ID)==T),which(colnames(A)%in%c('ID','logFC','Trend'))]
    BAC<-B[((B$ID%in%A$ID)==T & (B$ID%in%C$ID)==T),which(colnames(B)%in%c('ID','logFC','Trend'))]
    CAB<-C[((C$ID%in%A$ID)==T & (C$ID%in%B$ID)==T),which(colnames(C)%in%c('ID','logFC','Trend'))]
    ABC<-merge(ABC,BAC,by='ID')
    ABC<-merge(ABC,CAB,by='ID')
    rownames(ABC)<-ABC$ID
    ABC<-ABC[,-1]
  }
  A_only<-A[((A$ID%in%B$ID)==F & (A$ID%in%C$ID)==F),which(colnames(A)%in%c('ID','logFC','Trend'))]
  rownames(A_only)<-A_only$ID
  A_only<-A_only[,-1]
  B_only<-B[((B$ID%in%A$ID)==F & (B$ID%in%C$ID)==F),which(colnames(A)%in%c('ID','logFC','Trend'))]
  rownames(B_only)<-B_only$ID
  B_only<-B_only[,-1]
  C_only<-C[((C$ID%in%A$ID)==F & (C$ID%in%B$ID)==F),which(colnames(A)%in%c('ID','logFC','Trend'))]
  rownames(C_only)<-C_only$ID
  C_only<-C_only[,-1]
  UP<-c(UP,sum(A_only$Trend=='UP'));DOWN<-c(DOWN,sum(A_only$Trend=='DOWN'));Change<-c(Change,sum(A_only$Trend=='Change'))
  UP<-c(UP,sum(B_only$Trend=='UP'));DOWN<-c(DOWN,sum(B_only$Trend=='DOWN'));Change<-c(Change,sum(B_only$Trend=='Change'))
  UP<-c(UP,sum(C_only$Trend=='UP'));DOWN<-c(DOWN,sum(C_only$Trend=='DOWN'));Change<-c(Change,sum(C_only$Trend=='Change'))
  if (dim(AB)[1]==0){
    OvAB<-data.frame()
    UP<-c(UP,0);DOWN<-c(DOWN,0);Change<-c(Change,0)
  }else{
    tmp<-NULL
    for (t in 1:dim(AB)[1]) tmp<-c(tmp,ifelse(AB$Trend.x[t]==AB$Trend.y[t],AB$Trend.x[t],'Change'))
    OvAB<-data.frame(logFC_A=AB$logFC.x,logFC_B=AB$logFC.y,Trend=tmp)
    rownames(OvAB)<-rownames(AB)
    AB<-OvAB[order(OvAB$Trend),]
    UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
  }
  if (dim(AC)[1]==0){
    OvAc<-data.frame()
    UP<-c(UP,0);DOWN<-c(DOWN,0);Change<-c(Change,0)
  }else{
    tmp<-NULL
    for (t in 1:dim(AC)[1]) tmp<-c(tmp,ifelse(AC$Trend.x[t]==AC$Trend.y[t],AC$Trend.x[t],'Change'))
    OvAC<-data.frame(logFC_A=AC$logFC.x,logFC_C=AC$logFC.y,Trend=tmp)
    rownames(OvAC)<-rownames(AC)
    AC<-OvAC[order(OvAC$Trend),]
    UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
  }
  if (dim(BC)[1]==0){
    OvBC<-data.frame()
    UP<-c(UP,0);DOWN<-c(DOWN,0);Change<-c(Change,0)
  }else{
    tmp<-NULL
    for (t in 1:dim(BC)[1]) tmp<-c(tmp,ifelse(BC$Trend.x[t]==BC$Trend.y[t],BC$Trend.x[t],'Change'))
    OvBC<-data.frame(logFC_B=BC$logFC.x,logFC_C=BC$logFC.y,Trend=tmp)
    rownames(OvBC)<-rownames(BC)
    BC<-OvBC[order(OvBC$Trend),]
    UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
  }
  if (dim(ABC)[1]==0){
    OvABC<-data.frame()
    UP<-c(UP,0);DOWN<-c(DOWN,0);Change<-c(Change,0)
  }else{
    tmp<-NULL
    for (t in 1:dim(ABC)[1]) tmp<-c(tmp,ifelse(((ABC$Trend.x[t]==ABC$Trend.y[t]) & (ABC$Trend.x[t]==ABC$Trend[t])),ABC$Trend.x[t],'Change'))
    OvABC<-data.frame(logFC_A=ABC$logFC.x,logFC_B=ABC$logFC.y,logFC_C=ABC$logFC,Trend=tmp)
    rownames(OvABC)<-rownames(ABC)
    ABC<-OvABC[order(OvABC$Trend),]
    UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
  }
  counts<-data.frame(Contrast=c('A_only','B_only','C_only','AB','AC','BC','ABC'),Count=c(dim(A_only)[1],dim(B_only)[1],dim(C_only)[1],dim(AB)[1],dim(AC)[1],dim(BC)[1],dim(ABC)[1]),UP=UP,DOWN=DOWN,Change=Change)
  venn<-list(A_only=A_only,B_only=B_only,C_only=C_only,AB=AB,BC=BC,AC=AC,ABC=ABC)
  return(list(venn_df=counts,table=venn))
}

# Caluclate overlap for two lists
get_overlap2<-function(A,B){
  colnames(A)<-c('ID','logFC')
  colnames(B)<-c('ID','logFC')
  UP<-NULL;DOWN<-NULL;Change<-NULL
  if (class(A$logFC)!='numeric'){
    A$logFC<-gsub(",", ".", gsub("\\.", "", A$logFC))
    A$Trend<-sapply(as.numeric(A$logFC), function(x) ifelse(x > 0,'UP','DOWN')) 
  }else{ A$Trend<-sapply(A$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
  if (class(B$logFC)!='numeric'){
    B$logFC<-gsub(",", ".", gsub("\\.", "", B$logFC))
    B$Trend<-sapply(as.numeric(B$logFC), function(x) ifelse(x > 0,'UP','DOWN')) 
  }else{ B$Trend<-sapply(B$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
  AB<-A[(A$ID%in%B$ID)==T,which(colnames(A)%in%c('ID','logFC','Trend'))]
  BA<-B[(B$ID%in%A$ID)==T,which(colnames(B)%in%c('ID','logFC','Trend'))]
  A_only<-A[(A$ID%in%B$ID)==F,which(colnames(A)%in%c('ID','logFC','Trend'))]
  B_only<-B[(B$ID%in%A$ID)==F,which(colnames(B)%in%c('ID','logFC','Trend'))]
  AB<-merge(AB,BA,by='ID')
  UP<-c(UP,sum(A_only$Trend=='UP'));DOWN<-c(DOWN,sum(A_only$Trend=='DOWN'));Change<-c(Change,sum(A_only$Trend=='Change'))
  UP<-c(UP,sum(B_only$Trend=='UP'));DOWN<-c(DOWN,sum(B_only$Trend=='DOWN'));Change<-c(Change,sum(B_only$Trend=='Change'))
  rownames(A_only)<-A_only$ID
  A_only<-A_only[,-1]
  A_only<-A_only[order(A_only$Trend),]
  rownames(B_only)<-B_only$ID
  B_only<-B_only[,-1]
  B_only<-B_only[order(B_only$Trend),]
  tmp<-NULL
  for (t in 1:dim(AB)[1]) tmp<-c(tmp,ifelse(AB$Trend.x[t]==AB$Trend.y[t],AB$Trend.x[t],'Change'))
  OvAB<-data.frame(logFC_A=AB$logFC.x,logFC_B=AB$logFC.y,Trend=tmp)
  rownames(OvAB)<-AB$ID
  AB<-OvAB[order(OvAB$Trend),]
  UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
  counts<-data.frame(Contrast=c('A_only','B_only','AB'),Count=c(dim(A_only)[1],dim(B_only)[1],dim(AB)[1]),UP=UP,DOWN=DOWN,Change=Change)
  venn<-list(A_only=A_only,B_only=B_only,AB=AB)
  return(list(venn_df=counts,table=venn,dim=c(dim(A)[1],dim(B)[1])))
}