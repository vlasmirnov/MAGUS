'''
Created on May 14, 2020

@author: Vlad
'''

import os
import time


from align.merge.graph_builder import buildGraph
from align.merge.graph_cluster.clusterer import clusterGraph
from align.merge.graph_trace.tracer import findTrace
from align.merge.optimizer import optimizeTrace
from configuration import Configs

def mergeSubalignments(task):
    Configs.log("Merging {} subaligments..".format(len(task.subalignmentPaths)))
    time1 = time.time()  
        
    buildGraph(task)
    clusterGraph(task.graph)
    findTrace(task.graph)
    
    if max(task.graph.subalignmentLengths) > 10000:
        task.graph.cheaterAlignment(task.outputFile)
    else:
        task.graph.addSingletonClusters()
        optimizeTrace(task.graph)
        task.graph.clustersToAlignment(task.outputFile)
    
    time2 = time.time()  
    Configs.log("Merged {} subalignments into {} in {} sec..".format(len(task.subalignmentPaths), task.outputFile, time2-time1))
