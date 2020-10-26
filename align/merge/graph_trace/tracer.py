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
from align.merge.graph_trace.rg_search import rgSearch
from align.merge.graph_trace.rg_fast_search import rgFastSearch
from align.merge.graph_trace.naive import naiveClustering

def findTrace(graph):
    time1 = time.time() 
    
    if Configs.graphTraceMethod == "minclusters":
        graph.clusters = minClustersSearch(graph)       
    elif Configs.graphTraceMethod == "fm":            
        graph.clusters = fmAlgorithm(graph)        
    elif Configs.graphTraceMethod == "mwtgreedy":
        graph.clusters = mwtGreedySearch(graph)
    elif Configs.graphTraceMethod == "mwtsearch":
        graph.clusters = mwtSearch(graph)
    elif Configs.graphTraceMethod == "rg":
        graph.clusters = rgSearch(graph)
    elif Configs.graphTraceMethod == "rgfast":
        graph.clusters = rgFastSearch(graph)
    elif Configs.graphTraceMethod == "naive":
        graph.clusters = naiveClustering(graph)
    
    #for cluster in graph.clusters:
    #    Configs.log("{}".format(cluster.sorted()))
    
        
    time2 = time.time()
    Configs.log("Found alignment graph trace in {} sec..".format(time2-time1))
    Configs.log("Found a trace with {} clusters and a total cost of {}".format(len(graph.clusters), graph.computeClusteringCost(graph.clusters)))
    
    
    
    