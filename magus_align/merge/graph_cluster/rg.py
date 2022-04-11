'''
Created on Aug 23, 2020

@author: Vlad
'''

from magus_configuration import Configs
from magus_align.merge.graph_trace.rg_search import rgCluster
from magus_align.merge.graph_trace.rg_fast_search import rgFastCluster

def rgClustering(graph):
    Configs.log("Building a region-growing graph clustering..")
    
    k = len(graph.context.subalignments)
    lowerBound = [graph.subsetMatrixIdx[i] for i in range(k)]
    upperBound = [graph.subsetMatrixIdx[i] + graph.subalignmentLengths[i] for i in range(k)] 
    graph.clusters =  rgCluster(graph, lowerBound, upperBound, False)
    graph.writeClustersToFile(graph.clusterPath)

def rgFastClustering(graph):
    Configs.log("Building a fast region-growing graph clustering..")
    
    k = len(graph.context.subalignments)
    lowerBound = [graph.subsetMatrixIdx[i] for i in range(k)]
    upperBound = [graph.subsetMatrixIdx[i] + graph.subalignmentLengths[i] for i in range(k)] 
    graph.clusters =  rgFastCluster(graph, lowerBound, upperBound, False)
    graph.writeClustersToFile(graph.clusterPath)
