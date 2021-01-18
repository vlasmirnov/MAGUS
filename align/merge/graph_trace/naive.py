'''
Created on Aug 23, 2020

@author: Vlad
'''


from configuration import Configs

def atomizedClustering(graph):
    Configs.log("Building a fully atomized clustering..")
    
    k = len(graph.context.subalignments)
    lowerBound = [graph.subsetMatrixIdx[i] for i in range(k)]
    upperBound = [graph.subsetMatrixIdx[i] + graph.subalignmentLengths[i] for i in range(k)] 
    graph.clusters = atomizedCluster(lowerBound, upperBound)

def naiveClustering(graph):
    Configs.log("Building a naive left-justified clustering..")
    
    k = len(graph.context.subalignments)
    lowerBound = [graph.subsetMatrixIdx[i] for i in range(k)]
    upperBound = [graph.subsetMatrixIdx[i] + graph.subalignmentLengths[i] for i in range(k)] 
    graph.clusters = naiveCluster(lowerBound, upperBound)

def atomizedCluster(lowerBound, upperBound):
    clusters = []
    for j in range(len(lowerBound)):
        for i in range(lowerBound[j], upperBound[j]):
            clusters.append([i])
    return clusters

def naiveCluster(lowerBound, upperBound):
    clusters = []
    i = 0
    while True:
        cluster = []
        for j in range(len(lowerBound)):
            if lowerBound[j] + i < upperBound[j]:
                cluster.append(lowerBound[j] + i)
        if len(cluster) == 0:
            break   
        clusters.append(cluster)     
        i = i+1
    return clusters