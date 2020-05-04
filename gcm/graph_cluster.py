'''
Created on Apr 14, 2020

@author: Vlad
'''

import time
import os

from configuration import Configs
from tools import external_tools

def clusterGraph(graph):
    time1 = time.time()
    
    if not os.path.exists(Configs.clusterPath):
        external_tools.runMcl(Configs.graphPath, Configs.mclInflationFactor, Configs.workingDir, Configs.clusterPath)
    graph.readClustersFromFile(Configs.clusterPath)
    
    time2 = time.time()  
    Configs.log("Clustered the graph in {} sec..".format(time2-time1))
    
    purgeDuplicateClusters(graph)
    purgeClusterRedundancies(graph)
        
    time3 = time.time()  
    Configs.log("Purged the clusters in {} sec..".format(time3-time2))
    
def purgeDuplicateClusters(graph):
    uniqueClusters = set()
    newclusters = []
    for cluster in graph.clusters:
        cluster.sort()
        clusterTuple = tuple(cluster)
        if clusterTuple not in uniqueClusters:
            uniqueClusters.add(clusterTuple)
            newclusters.append(cluster)
    graph.clusters = newclusters
    Configs.log("Purged duplicate clusters. Found {} unique clusters..".format(len(graph.clusters)))

def purgeClusterRedundancies(graph):
    redundantCols = {}
    redundantRows = {}
    elementScores = {}
    for a, cluster in enumerate(graph.clusters):
        for b in cluster:
            bsub, bpos = graph.matSubPosMap[b] 
            redundantCols[a, bsub] = redundantCols.get((a, bsub), []) + [(a, b)] 
            redundantRows[b] = redundantRows.get(b, []) + [(a, b)]
            
            scoresum = 0
            for c in cluster:
                csub, cpos = graph.matSubPosMap[c]
                if bsub != csub:
                    scoresum  = scoresum  + graph.matrix[b].get(c,0)
            elementScores[a, b] = scoresum 
    
    problemCols = [(a,b) for a,b in redundantCols if len(redundantCols[a,b]) > 1]
    problemRows = [a for a in redundantRows if len(redundantRows[a]) > 1]
    Configs.log("Found {} row redundancies and {} column redundancies..".format(len(problemRows), len(problemCols)))
    
    sortedScores = list(elementScores.keys())
    sortedScores.sort(key = lambda x : elementScores[x])
    
    for a,b in sortedScores:
        bsub, bpos = graph.matSubPosMap[b] 
        if len(redundantCols[a, bsub]) > 1 or len(redundantRows[b]) > 1:
            graph.clusters[a].remove(b)
            redundantCols[a, bsub].remove((a,b))
            redundantRows[b].remove((a,b))                
    
    problemCols = [(a,b) for a,b in redundantCols if len(redundantCols[a,b]) > 1]
    problemRows = [a for a in redundantRows if len(redundantRows[a]) > 1]
    Configs.log("Finished redundancy sweep. Now {} row redundancies and {} column redundancies..".format(len(problemRows), len(problemCols)))
    
    graph.clusters = [cluster for cluster in graph.clusters if len(cluster) > 1]
    Configs.log("Purged row/column redundancies. Found {} proper clusters..".format(len(graph.clusters)))