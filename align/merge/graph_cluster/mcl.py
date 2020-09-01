'''
Created on Apr 14, 2020

@author: Vlad
'''

import time
import os

from configuration import Configs
from tools import external_tools


def runMclClustering(graph):  
    Configs.log("Running MCL alignment graph clustering..")
    clusterPath = os.path.join(graph.workingDir, "clusters.txt")
    
    if not os.path.exists(clusterPath):
        external_tools.runMcl(graph.graphPath, Configs.mclInflationFactor, graph.workingDir, clusterPath)
    else:
        Configs.log("Found existing cluster file {}".format(clusterPath))
    graph.readClustersFromFile(clusterPath)
    
    
