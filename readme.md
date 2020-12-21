# Tic2D

## Overview

Tic2D can be used for cell clustering and pseudo-temporal reconstruction of RNA-seq datasets. Compared with other trajectory inference algorithms, Tic2D has higher general applicability and prediction accuracy in determining cell state and inferring cell trajectory. In addition, the reconstructed trajectory allows us to identify the key genes involved in cell biology and understand their expression patterns and functional roles in the dynamic progression of cells.
Tic2D uses k-nearest neighbor graph and louvain algorithm to group genes and cells to construct a binary matrix, and then perform clustering on the binary matrix to obtain cluster labels.
Tic2D constructs an MST to connect all cluster centers, and uses the Principal curve to obtain a smooth trajectory as the cell development time axis, and then projects the cells to this axis to sort the cells.

### All the data and code need to be downloaded to run ###

### additional is some Python version code we used, Here is for reference. ###

### Please install dynutils, igraph, princurve, Rtsne before use. ###

 ### re_matrix()  and Trajectory_Infer()  are functions in the R file. ###

re_matrix(data,k) : The clustering method used in our experiments for R.

#------EXAMPLE------#

#k:k is the number of clusters
#Binarymatrix:Binarymatrix is the generated binary matrix
#which can be used for trajectory inference

exprdata <- "HSMM_0.15.txt"
procdata <- read.table(exprdata, sep='\t', header=T,row.names = 1, quote = "", comment="", check.names=F)
data <- as.matrix(procdata)
k <- 4
Binarymatrix<-re_matrix(data,k)



Trajectory_Infer(X,kmeans_clust,k,startCluster,exprsub): Get the cell trajectories and pseudo-time of cells. 

#k:k is the number of clusters
#startCluster:startCluster is the starting cluster of the inferred trajectory

#------EXAMPLE------#

X<-Rtsne(Binarymatrix,dims = 2, pca = FALSE,perplexity = 30,theta = 0.0,max_iter = 1000, check_duplicates = FALSE)
X<-X$Y
kmeans_clust <- stats::kmeans(as.matrix(X),centers = k)
k<-7
startCluster<-1
result<-Trajectory_Infer(X,kmeans_clust,k,startCluster,exprsub)