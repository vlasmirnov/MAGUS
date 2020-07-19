'''
Created on Jul 17, 2020

@author: Vlad
'''

import os
import time

from configuration import Configs
from align.merge.graph_trace.min_clusters import minClustersSearch
from align.merge.graph_trace.fm import fmAlgorithm
from align.merge.graph_trace.mwt_search import mwtGreedySearch, mwtSearch

def findTrace(graph):
    time1 = time.time() 
    
    if Configs.graphTraceMethod == "minclusters":
        graph.clusters = minClustersSearch(graph)
    
    else:    
        if Configs.graphClusterMethod == "mcl":
            graph.buildNodeEdgeDataStructureFromClusters()
        else:
            graph.buildNodeEdgeDataStructure()
        
        if Configs.graphTraceMethod == "fm":
            graph.clusters = fmAlgorithm(graph)
        
        elif Configs.graphTraceMethod == "mwtgreedy":
            graph.clusters = mwtGreedySearch(graph)
        
        elif Configs.graphTraceMethod == "mwtsearch":
            graph.clusters = mwtSearch(graph)
        
    time2 = time.time()
    Configs.log("Found alignment graph trace in {} sec..".format(time2-time1))
    Configs.log("Found a trace with {} clusters and a total cost of {}..".format(len(graph.clusters), graph.computeClusteringCost(graph.clusters)))