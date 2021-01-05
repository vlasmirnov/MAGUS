'''
Created on Aug 23, 2020

@author: Vlad
'''

import os
import time

from configuration import Configs

from align.merge.graph_cluster.mcl import runMclClustering
from align.merge.graph_cluster.mlr_mcl import runMlrMclClustering
from align.merge.graph_cluster.rg import rgClustering


def clusterGraph(graph):
    time1 = time.time()
    
    if os.path.exists(graph.clusterPath):
        Configs.log("Found existing cluster file {}".format(graph.clusterPath))
        graph.readClustersFromFile(graph.clusterPath)
        
    elif Configs.graphClusterMethod == "mcl":
        runMclClustering(graph)      
        
    elif Configs.graphClusterMethod == "mlrmcl":
        runMlrMclClustering(graph)
        
    elif Configs.graphClusterMethod == "rg":
        rgClustering(graph)
        
    else:
        Configs.log("No alignment graph clustering requested..")
    
    time2 = time.time()  
    Configs.log("Clustered the graph in {} sec..".format(time2-time1))
