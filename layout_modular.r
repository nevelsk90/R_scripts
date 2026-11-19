### a custom layout function to group nodes by communities 
### modified from https://stackoverflow.com/questions/16390221/how-to-make-grouped-layout-in-igraph



layout.modular <- function(G, c , sel_com)
  {
# G - your graph
# c-  community information for every node
# sel_com - IDs of selected communtities 

  nm <- length(levels(as.factor(c$membership)))
  gr <- 2
  while(gr^2<nm){
    gr <- gr+1
  }
  i <- j <- 0
  for(cc in levels(as.factor(c$membership))){
    F <- delete.vertices(G,c$membership!=cc)
    F$layout <- layout_in_circle(F)
    # scale circles if they are selected or not
    if(cc%in%sel_com) {
      scSize <- sqrt(gorder(F))*0.1
      F$layout <- layout.norm(F$layout, i, i+scSize, j, j+scSize )
    } 
    else  {
      scSize <- sqrt(gorder(F))*0.05
      
      F$layout <- layout.norm(F$layout, i , i+scSize, j, j+scSize )
    }
    
    G$layout[c$membership==cc,] <- F$layout
    if(i>=gr){
      i <- 0
      if(j>=gr){
        j <- 0
      }else{
        j <- j+1.5
      }
    } else{
      i <- i+scSize*1.4
    }
  }
  return(G$layout)
}

#### how to use
# find communities
	net.com.edge <- cluster_edge_betweenness( graph_toPlot , modularity=T, directed = F, weights = NULL)
# select only communities with e.g. more than 30 members  
 	sel_com <- which(summary(as.factor( net.com.edge$membership))>30)
# create a custom layout to group nodes by communities on a plot! 
	graph_toPlot$layout <- layout_with_mds(graph_toPlot)
#  use the custom layout function to group nodes by communities 
	l <- layout.modular( G=graph_toPlot , c=net.com.edge , sel_com=sel_com) 
# make a plot
	plot( graph_toPlot , layout=l )
