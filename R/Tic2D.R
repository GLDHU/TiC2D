library(TSCAN)
library(dynutils)
library(igraph)
library(princurve)
library(Rtsne)
library(flexclust)

#example
#exprdata <- "data.txt"
#exprsub<-"real_cluster.txt"
#procdata <- read.table(exprdata, sep='\t', header=T,row.names = 1, quote = "", comment="", check.names=F)
#data <- as.matrix(procdata)
#k <- 4
#Binarymatrix<-re_matrix(data,k)

#X<-Rtsne(Binarymatrix,dims = 2, pca = FALSE,perplexity = 30,theta = 0.0,max_iter = 1000, check_duplicates = FALSE)
#X<-X$Y
#kmeans_clust <- stats::kmeans(as.matrix(X),centers = k)
#k<-4
#startCluster<-2
#result<-Trajectory_Infer(X,kmeans_clust,k,startCluster)


re_matrix<- function(data, k){
  Grouping<-function(data, k){
    distance_matrix <- as.matrix(dist(data))
    distance_matrix1 <- data.frame(distance_matrix)
    knn_ids <- distance_matrix
    for(i in seq_len(length(distance_matrix1))){
      knn_ids[i,] <- order(distance_matrix[i,])
    }
    aa <- c()
    for(i in seq_len(length(distance_matrix1))){
      k1 <- 0
      for(jj in seq_len(length(distance_matrix1))){
        j <- knn_ids[i, jj]
        if(i != j){
          aa<- c(aa, i, j)
          k1 <- k1+1
        }
        if (k1 > k){
          break
        }
      }
    }
    g <- make_empty_graph(n = length(distance_matrix1), directed = FALSE) %>%
      add_edges(aa) %>%
      set_edge_attr("color", value = "red") 
    #plot(g)
    fen <- cluster_louvain(g)
    return(fen)
  }
  cluster_gene<-Grouping(data,k)
  d<-membership(cluster_gene)
  d1<-as.vector(d)
  data2<-cbind(data,d1)
  max_n<-max(d1)
  ncell<-ncol(data2)
  cell_matrix<-matrix(0,nrow=ncell-1)
  for(ii in seq_len(max_n)){
    temp1<-data2[data2[,ncol(data2)]==ii,1:(ncol(data2)-1)]
    temp2<-t(temp1)
    cluster_cell<-Grouping(temp2,k)
    b<-membership(cluster_cell)
    b1<-as.vector(b)
    max_b1<-max(b1)
    ncell<-nrow(temp2)
    cell_matrix_temp<-matrix(0,nrow=ncell,ncol=max_b1)
    for(j in seq_len(ncell)){
      cell_matrix_temp[j,b1[j]]=1
    }
    cell_matrix<-cbind(cell_matrix,cell_matrix_temp)
  }
  cell_matrix<-as.matrix(cell_matrix[,-1])
  Binarymatrix<-cell_matrix
  Binarymatrix<-as.matrix(Binarymatrix)
  X<-Rtsne(Binarymatrix,dims = 2, pca = FALSE,perplexity = 30,theta = 0.0,max_iter = 1000, check_duplicates = FALSE)
  X<-X$Y
  kmeans_clust <- stats::kmeans(as.matrix(X),centers = k)
  labels_pred<-kmeans_clust$cluster
  write.csv(labels_pred,"labels_pred.csv")
  #result<-Trajectory_Infer(X,kmeans_clust,k,startCluster)
  return(Binarymatrix)
}
Trajectory_Infer<-function(X,kmeans_clust,k,startCluster){
  #set.seed(12345)
  centers <- kmeans_clust$centers
  clusterLabels <- as.matrix(kmeans_clust$cluster)
  X.original <- X
  clusters <- c(1:k)
  nclus <- k
  start.clus <- as.character(startCluster)
  dist.fun<-NULL
  # calculate the euclidean space between clusters
  eucl_dist <- as.matrix(stats::dist(centers))
  knn_distances <- function(dist, k, self_loops=FALSE) {
    knn(dist, k, self_loops = self_loops)$distances
  }
  knn <- function(dist, k, self_loops = FALSE) {
    # input checks  
    if (class(dist) == "dist")
      dist <- as.matrix(dist)
    if (nrow(dist) < 2)
      stop(sQuote("dist"), " needs to consist of at least 2 rows")
    # initialise matrices with NAs   
    indices <- knndist <- matrix(
      NA,
      nrow = nrow(dist),
      ncol = k,
      dimnames = list(
        rownames(dist),
        if (k == 0) c() else paste0("knn", seq_len(k))
      )
    )
    # use dimnames if possible
    row_ix <- if (!is.null(rownames(dist))) rownames(dist) else seq_len(nrow(dist))
    col_ix <- if (!is.null(colnames(dist))) colnames(dist) else seq_len(ncol(dist))
    # fill matrix by sample
    if (self_loops) {
      for (i in row_ix) {
        indices[i,] <- utils::head(order(dist[i,]), k)
        knndist[i,] <- dist[i,indices[i,]]
      }
    } else {
      diag(dist) <- 0 # just to make sure
      for (i in row_ix) {
        indices[i,] <- head(order(dist[i,]), k+1)[-1]
        knndist[i,] <- dist[i,indices[i,]]
      }
    }         
    # return KNN distances
    list(
      indices = indices,
      distances = knndist
    )
  }
  # cal''culate the densities along the straight lines between any two cluster centers
  density_dist <- sapply(seq_len(k), function(i) {
    sapply(seq_len(k), function(j) {
      if (i == j) {
        0
      } else {
        twocent <- centers[c(i,j), , drop = FALSE]
        segment_pts <- apply(twocent, 2, function(x) seq(x[[1]], x[[2]], length.out = 20))
        dists <- euclidean_distance(segment_pts, X)
        mean(knn_distances(dists, 10, self_loops=TRUE))
      }
    })
  })
  # combine both distance matrices  
  D <- eucl_dist * density_dist
  # Minimum spanning tree constructed for identifing lineages
  mstree <- ape::mst(D)
  forest <- mstree
  lineages <- list()
  # identify sub-trees 
  subtrees <- subtrees.update <- forest
  diag(subtrees) <- 1
  while(sum(subtrees.update) > 0){
    subtrees.new <- apply(subtrees,2,function(col){
      rowSums(subtrees[,as.logical(col), drop=FALSE]) > 0
    })
    subtrees.update <- subtrees.new - subtrees
    subtrees <- subtrees.new
  }
  subtrees <- unique(subtrees)
  trees <- lapply(seq_len(nrow(subtrees)),function(ri){
    colnames(forest)[subtrees[ri,]]
  })
  trees <- trees[order(vapply(trees,length,0),decreasing = TRUE)]
  ntree <- length(trees)
  # identify lineages (paths through trees)
  for(tree in trees){
    if(length(tree) == 1){
      lineages[[length(lineages)+1]] <- tree
      next
    }
    tree.ind <- rownames(forest) %in% tree
    tree.graph <- forest[tree.ind, tree.ind, drop = FALSE]
    degree <- rowSums(tree.graph)
    g <- graph.adjacency(tree.graph, mode="undirected")
    # if you have starting cluster(s) in this tree, draw lineages
    # to each leaf
    if(sum(start.clus %in% tree) > 0){
      starts <- start.clus[start.clus %in% tree]
      ends <- rownames(tree.graph)[
        degree == 1 & ! rownames(tree.graph) %in% starts]
      for(st in starts){
        paths <- shortest_paths(g, from = st, to = ends, 
                                mode = 'out', 
                                output = 'vpath')$vpath
        for(p in paths){
          lineages[[length(lineages)+1]] <- names(p)
        }
      }
    }
  }
  #Visualization of the minimum spanning tree
  plot(X[,1], X[,2], col = kmeans_clust$cluster, xlab ="", ylab ="")
  #points(centers, col = 1:k, pch = 16, cex=2)
  #text(centers,labels=1:nclus,pos="1",offset=-1,col="blue")
  dims = seq_len(2)
  for(i in seq_len(nclus-1)){
    for(j in seq(i+1,nclus)){
      if(mstree[i,j] == 1){
        lines(centers[c(i,j), dims], lwd = 2, col = "black")
      }
    }
  }
  plot(X[,1], X[,2], col = kmeans_clust$cluster, xlab ="", ylab ="")
  #Visualization of the trajectory curves
  Trajectory <- NULL
  X1 <- cbind.data.frame(X,clusterLabels)
  orders<-list()
  #View(lineages)
  for(i in seq_len(length(lineages))){
    c = data.frame(v1 = 0,v2 = 0)
    c = c[-1,]
    for(j in seq(length(lineages[[i]]))){
      nu <- as.integer(lineages[[i]][j])
      c1 <- X1[clusterLabels==nu,]
      c <- rbind(c, c1[,-3])
    }
    c <- as.matrix(c) 
    fit <- principal_curve(c)
    s = fit$s[fit$ord, ,drop = FALSE]
    fit <- project_to_curve(X, s)
    lines(fit, lwd = 2, col = "black")
    Trajectory[[i]] <- fit
    
    #order<-fit$lambda
    #cell<-subpopulation$cell
    #cellorder<-cbind(data.frame(cell), data.frame(order))
    #A<-cellorder[order(cellorder$order,decreasing=F),]
    #order_1<-A$cell
    #orders[i]<-list(order_1)
  }
  #print(orderscore(subpopulation, orders))
  return(Trajectory)
}

