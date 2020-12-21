# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from sklearn import metrics
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
import networkx as nx
from sklearn.metrics.pairwise import pairwise_distances
import community
from sklearn.cluster import KMeans
from time import *
from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.decomposition import PCA

#Load binary matrix
def _loaddata_matrix(fname):
    data_matrix_original = []
    instance_names = []
    gene_names = []
    with open(fname) as f:
        for i, line in enumerate(f):
            if i == 0:
                instance_names = line.strip().split()[0:]
            if i > 0:
                tokens = line.strip().split('\t')
                gene_names.append(tokens[0])
                value_list = tokens[1:]
                vals = [float(j) for j in value_list]
                data_matrix_original.append(vals)
    data_matrix = np.array(data_matrix_original).T
    rows, cols = data_matrix.shape
    return data_matrix, gene_names, instance_names

#Load transposed binary matrix
def _loaddata_matrix1(fname):
    data_matrix_original = []
    instance_names = []
    gene_names = []
    with open(fname) as f:
        for i, line in enumerate(f):
            if i == 0:
                instance_names = line.strip().split()[0:]
            if i > 0:
                tokens = line.strip().split('\t')
                gene_names.append(tokens[0])
                value_list = tokens[1:]
                vals = [float(j) for j in value_list]
                data_matrix_original.append(vals)
    data_matrix = np.array(data_matrix_original)
    rows, cols = data_matrix.shape
    return data_matrix, gene_names, instance_names


def loadDataSet(fileName, curLine=None):  #
    # type: (object) -> object
    dataMat = []
    data = []
    fr = open(fileName)

    for line in fr.readlines():
        curLine = line.strip("\n\r").split("\t")

        dataMat.append(curLine)
    dataMat = np.array(dataMat)
    #   dataMat = np.delete(dataMat, 0, 1)
    #  dataMat = np.delete(dataMat, 0, 0)
    for i in dataMat:
        b = [float(j) for j in i]
        data.append(b)
    dataMat = np.array(data)
    return dataMat


def loadDataSet2(fileName, curLine=None):  #
    # type: (object) -> object
    dataMat = []
    data = []
    fr = open(fileName)

    for line in fr.readlines():
        curLine = line.strip("\n\r").split("\t")

        dataMat.append(curLine)
    dataMat = np.array(dataMat)
    dataMat = np.delete(dataMat, 0, 1)
    dataMat = np.delete(dataMat, 0, 0)
    for i in dataMat:
        b = [float(j) for j in i]
        data.append(b)
    dataMat = np.array(data).T
    return dataMat


def loadDataSet1(fileName, curLine=None):  #
    # type: (object) -> object
    dataMat = []
    data = []
    fr = open(fileName)

    for line in fr.readlines():
        curLine = line.strip("\n\r").split("\t")

        dataMat.append(curLine)
    dataMat = np.array(dataMat)
    # dataMat = np.delete(dataMat, 0, 1)
    # dataMat = np.delete(dataMat, 0, 0)
    for i in dataMat:
        b = [float(j) for j in i]
        data.append(b)
    dataMat = np.array(data)
    return dataMat

#Set knn graph
def _add_knn_links(self, graph, target, nneighbors_th=1, kernel_matrix=None, knn_ids=None):
    size = kernel_matrix.shape[0]
    for i in range(size):
        # add edges to the k nns with same class
        k = 0
        for jj in range(size):
            j = knn_ids[i, jj]
            if i != j:
                if target[i] == target[j]:
                    graph.add_edge(i, j, edge_type='knn', rank=k)
                    k += 1
            # after having added at most nneighbors_th links exit
            if k > nneighbors_th:
                break
    return graph


# k k is the number of genes
# n is the number of rows
def Binarymatrix(k, labels, n):
    Binmat = np.mat(np.zeros((n, k)))
    label = labels
    for i in range(n):
        for j in range(k):
            if (int(label[i]) == j):
                Binmat[i, j] = 1
                break
    return Binmat



def labellen(old_list):
    newList = []
    for x in old_list:
        if x not in newList:
            newList.append(x)
    return len(newList)


# knn-graph + louvain to cluster genes, get n subdatasets, 'k1' is k in knn
#return name is the number of subsets name1 is the number of cells
def GenesCluster(data, k1):
    data_matrix, gene_names, instance_names = _loaddata_matrix1(data)
    gene_names = np.array(gene_names)
    instance_names = np.array(instance_names)
    size = data_matrix.shape[0]
    distance_matrix = pairwise_distances(data_matrix)
    knn_ids = np.argsort(distance_matrix)
    nneighbors_th = k1
    graph = nx.Graph()
    for v in range(size):
        graph.add_node(v, outlier=False)
    for i in range(size):
        # add edges to the k nns with same class
        k = 0
        for jj in range(size):
            j = knn_ids[i, jj]
            if i != j:
                graph.add_edge(i, j, edge_type='knn', rank=k)
                k += 1
            # after having added at most nneighbors_th links exit
            if k > nneighbors_th:
                break
    for i in range(size):
        counter = 0
        # add edges to the knns
        for jj in range(1, int(nneighbors_th) + 1):
            j = knn_ids[i, jj]
            if i != j:
                # check that within the k-nn also i is a knn of j
                # i.e. use the symmetric nneighbor notion
                upto = int(nneighbors_th) + 1
                i_knns = knn_ids[j, :upto]
                if i in list(i_knns):
                    counter += 1
                    break
        if counter > 0:
            outlier_status = False
        else:
            outlier_status = True

        graph.node[i]['outlier'] = outlier_status
        # nx.draw(graph,node_size=30,with_labels = False)
         # plt.show()
    # nx.draw(graph, node_size=15, with_labels=False)
    # # nx.draw(graph, node_size=30, with_labels=False)
    # plt.title("the first figure", fontsize=20)
    # plt.show()
    partition = community.best_partition(graph,resolution=1.5)
    #drawing
    # size = float(len(set(partition.values())))
    # pos = nx.spring_layout(graph)
    # count = 0.
    # for com in set(partition.values()):
    #     count = count + 1.
    #     list_nodes = [nodes for nodes in partition.keys()
    #                   if partition[nodes] == com]
    #         nx.draw_networkx_nodes(graph, pos, list_nodes, node_size = 20,
    #                                    node_color = str(count / size))
    A = []
    for i in partition:
        A.append(partition[i])
    A = np.array(A)
    # nx.draw_networkx_edges(graph, pos, alpha=0.5)
    # plt.title("the second figure", fontsize=20)
    # plt.show()
    #print(metrics.silhouette_score(data_matrix, A, metric='euclidean'))
    C = A
    name1 = data_matrix.shape[1]
    #print(data_matrix.shape)
    a = np.c_[data_matrix, C]
    data2 = []
    data2 = np.append(instance_names, "cluster")
    dataExpr1 = pd.DataFrame(data=a, columns=data2, index=gene_names)
    #print(dataExpr1.shape)
    dataExpr1.to_csv('Genes classification for data.txt', sep='\t', encoding="utf-8", index=True, header=True)
    for name, group in dataExpr1.groupby("cluster"):
        a = group.drop(['cluster'], axis=1)
        a.to_csv('Genes classification for data genegroup' + str(name) + '.txt', sep='\t', encoding="utf-8", index=True,
                 header=True)
    return name, name1


# knn-graph + louvain to cluster cells for each subdatasets, get n binarymatrix
# 'n' is the number of subdatasets, 'k2' is k in knn  'c' is the number of cells
def Nbinarymatrix(n, k2, c):
    w = 0.0
    while w <= n:
        data_matrix, gene_names, instance_names = _loaddata_matrix(
            'Genes classification for data genegroup' + str(w) + '.txt')
        #print(instance_names)
        size = data_matrix.shape[0]
        distance_matrix = pairwise_distances(data_matrix)
        knn_ids = np.argsort(distance_matrix)
        nneighbors_th = k2
        graph = nx.Graph()
        for v in range(size):
            graph.add_node(v, outlier=False)
        for i in range(size):
            # add edges to the k nns with same class
            k = 0
            for jj in range(size):
                j = knn_ids[i, jj]
                if i != j:
                    graph.add_edge(i, j, edge_type='knn', rank=k)
                    k += 1
                # after having added at most nneighbors_th links exit
                if k > nneighbors_th:
                    break
        for i in range(size):
            counter = 0
            # add edges to the knns
            for jj in range(1, int(nneighbors_th) + 1):
                j = knn_ids[i, jj]
                if i != j:
                    # check that within the k-nn also i is a knn of j
                    # i.e. use the symmetric nneighbor notion
                    upto = int(nneighbors_th) + 1
                    i_knns = knn_ids[j, :upto]
                    if i in list(i_knns):
                        counter += 1
                        break
            if counter > 0:
                outlier_status = False
            else:
                outlier_status = True
            graph.node[i]['outlier'] = outlier_status
        #nx.draw(graph,node_size=30,with_labels = False)
        #plt.show()
        # nx.draw(graph, node_size=15, with_labels=False)
        # # nx.draw(graph, node_size=30, with_labels=False)
        # plt.title("the three figure", fontsize=20)
        # plt.show()
        partition = community.best_partition(graph,resolution=1.5)
        #drawing
        # size = float(len(set(partition.values())))
        # pos = nx.spring_layout(graph)
        # count = 0.
        # for com in set(partition.values()):
        #     count = count + 1.
        #     list_nodes = [nodes for nodes in partition.keys()
        #                   if partition[nodes] == com]
        #         nx.draw_networkx_nodes(graph, pos, list_nodes, node_size = 20,
        #                                    node_color = str(count / size))

        A = []
        for i in partition:
            A.append(partition[i])
        A = np.array(A)
        #       nx.draw_networkx_edges(graph, pos, alpha=0.5)
        #       plt.show()
        # nx.draw_networkx_edges(graph, pos, alpha=0.5)
        # plt.title("the four figure", fontsize=20)
        # plt.show()
        labelnum = int(labellen(A))
        smallmatrix = Binarymatrix(labelnum, A, c)
        np.savetxt('Genes classification for data genegroup' + str(w) + 'binarymatrix.txt', smallmatrix, fmt='%s',
                   delimiter='\t')
        w = w + 1

# merge n binarymatrix into a big Binarymatrix, perform Kmeans on Binarymatrix
def binarymatrixMerge(n, k):
    i = 0.0
    lmat = loadDataSet('Genes classification for data genegroup' + str(i) + 'binarymatrix.txt')
    #print(lmat.shape)
    i=1.0
    while i <= n:
        data = loadDataSet('Genes classification for data genegroup' + str(i) + 'binarymatrix.txt')
        #print(data.shape)
        lmat = np.hstack((lmat, data))
        i = i + 1
    #print(lmat.shape)
    np.savetxt('Binarymatrix.txt', lmat, fmt='%s', delimiter='\t')
    kmeans = KMeans(n_clusters=k, random_state=0).fit(lmat)
    labels_pred = kmeans.labels_
    # label = []
    # for i in labels_pred:
    #     if i == 0:
    #         i=k
    #     label.append(i)
    # #print(label)
    # label = np.array(label)
    return labels_pred



if __name__ == '__main__':
    file_data = "data.txt"
    file_time_true = "real_cluster.txt"
    k = 4
    name,name1 = GenesCluster(file_data, k)
    Nbinarymatrix(name, k, name1)
    #return labels_pred
    labels_pred = binarymatrixMerge(name, k)


    # labels_true = loadDataSet2(file_time_true)
    # labels_true = labels_true[0, :]
    # score = metrics.adjusted_rand_score(labels_true, labels_pred)
    # print(score)



