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
from align.merge.graph_cluster.clean_clusters import purgeClusterViolations, purgeDuplicateClusters

def clusterGraph(graph):
    time1 = time.time()
    
    if Configs.graphClusterMethod == "mcl": #or Configs.graphTraceMethod == "minclusters":
        runMclClustering(graph)
        #bins = max([len(c) for c in graph.clusters]) + 1
        #graph.plotClusterHistogram(graph.clusters, 70, os.path.join(graph.workingDir, "histogram.pdf"))
        purgeDuplicateClusters(graph)
        purgeClusterViolations(graph)
    elif Configs.graphClusterMethod == "mlrmcl":
        runMlrMclClustering(graph)
        #bins = max([len(c) for c in graph.clusters]) + 1
        #graph.plotClusterHistogram(graph.clusters, 70, os.path.join(graph.workingDir, "histogram.pdf"))
        purgeDuplicateClusters(graph)
        purgeClusterViolations(graph)
    elif Configs.graphClusterMethod == "rg":
        rgClustering(graph)
    else:
        Configs.log("No alignment graph clustering requested..")
    
    time2 = time.time()  
    Configs.log("Clustered the graph in {} sec..".format(time2-time1))
