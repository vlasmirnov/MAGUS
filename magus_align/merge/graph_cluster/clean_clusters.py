'''
Created on Aug 23, 2020

@author: Vlad
'''

from magus_configuration import Configs

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

def purgeClusterViolations(graph):
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
    Configs.log("Found {} row violations and {} column violations..".format(len(problemRows), len(problemCols)))
    
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
    Configs.log("Finished violations sweep. Now {} row violations and {} column violations..".format(len(problemRows), len(problemCols)))
    
    graph.clusters = [cluster for cluster in graph.clusters if len(cluster) > 1]
    Configs.log("Purged cluster violations. Found {} clean clusters..".format(len(graph.clusters)))
    