'''
Created on Jul 17, 2020

@author: Vlad
'''

import os
import time

from configuration import Configs
from align.merge.graph_cluster.clean_clusters import purgeClusterViolations, purgeDuplicateClusters
from align.merge.graph_trace.min_clusters import minClustersSearch
from align.merge.graph_trace.fm import fmAlgorithm
from align.merge.graph_trace.mwt_search import mwtGreedySearch, mwtSearch
from align.merge.graph_trace.rg_search import rgSearch
from align.merge.graph_trace.rg_fast_search import rgFastSearch
from align.merge.graph_trace.naive import naiveClustering

'''
Graph clusters must be refined into a "trace", a constrained clustering that corresponds to a valid MSA.
There are a variety of ways to do this, "minclusters" is usually the most dependable option.
"mwtgreedy" or "rgfast" might be used if there are scalability issues.
Some of these options don't require an existing clustering, and can work on the raw graph.
'''

def findTrace(graph):
    time1 = time.time() 
    
    if os.path.exists(graph.tracePath):
        Configs.log("Found existing trace file {}".format(graph.tracePath))
        graph.readClustersFromFile(graph.tracePath)
        
    else:
        purgeDuplicateClusters(graph)
        purgeClusterViolations(graph)
        
        if Configs.graphTraceMethod == "minclusters":
            minClustersSearch(graph)       
        elif Configs.graphTraceMethod == "fm":            
            fmAlgorithm(graph)        
        elif Configs.graphTraceMethod == "mwtgreedy":
            mwtGreedySearch(graph)
        elif Configs.graphTraceMethod == "mwtsearch":
            mwtSearch(graph)
        elif Configs.graphTraceMethod == "rg":
            rgSearch(graph)
        elif Configs.graphTraceMethod == "rgfast":
            rgFastSearch(graph)
        elif Configs.graphTraceMethod == "naive":
            naiveClustering(graph)
        
        graph.writeClustersToFile(graph.tracePath)
    
    
    time2 = time.time()
    Configs.log("Found alignment graph trace in {} sec..".format(time2-time1))
    Configs.log("Found a trace with {} clusters and a total cost of {}".format(len(graph.clusters), graph.computeClusteringCost(graph.clusters)))
    
    
    
    