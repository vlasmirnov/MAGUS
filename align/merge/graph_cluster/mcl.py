'''
Created on Apr 14, 2020

@author: Vlad
'''

import os

from configuration import Configs
from tools import external_tools
from helpers import tasks


def runMclClustering(graph):  
    Configs.log("Running MCL alignment graph clustering..")
    clusterPath = os.path.join(graph.workingDir, "clusters.txt")
    
    if not os.path.exists(clusterPath):
        task = external_tools.runMcl(graph.graphPath, Configs.mclInflationFactor, graph.workingDir, clusterPath)
        task.submitTask()
        task.waitForTask()
    else:
        Configs.log("Found existing cluster file {}".format(clusterPath))
    graph.readClustersFromFile(clusterPath)
    
    
