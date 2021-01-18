'''
Created on Apr 14, 2020

@author: Vlad
'''

from configuration import Configs
from tools import external_tools


def runMclClustering(graph):  
    Configs.log("Running MCL alignment graph clustering..")
    external_tools.runMcl(graph.graphPath, Configs.mclInflationFactor, graph.workingDir, graph.clusterPath).run()
    graph.readClustersFromFile(graph.clusterPath)
    
    
