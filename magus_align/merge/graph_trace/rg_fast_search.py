'''
Created on Aug 23, 2020

@author: Vlad
'''

import heapq
from collections import deque 

from magus_configuration import Configs

def rgFastSearch(graph):
    Configs.log("Finding graph trace with fast region-growing search..")
    
    k = len(graph.context.subalignments)
    lowerBound = [graph.subsetMatrixIdx[i] for i in range(k)]
    upperBound = [graph.subsetMatrixIdx[i] + graph.subalignmentLengths[i] for i in range(k)] 
    cuts = rgFastCluster(graph, lowerBound, upperBound, True)
    graph.clusters = cutsToClusters(graph, cuts)

def rgFastCluster(graph, lowerBound, upperBound, enforceTrace = True):
    initialCuts = initialSplit(graph, lowerBound, upperBound, enforceTrace)
    if len(initialCuts) == 2:
        return initialCuts
    #Configs.log("Starting with {} coarse cuts..".format(len(initialCuts)))
    cuts = []
    for i in range(len(initialCuts)-1):
        intervalCuts = rgFastCluster(graph, initialCuts[i], initialCuts[i+1], enforceTrace)
        cuts.extend(intervalCuts[:-1])
    
    cuts.append(list(upperBound))
    #Configs.log("Returning {} fine cuts..".format(len(cuts)))
    return cuts
        
def initialSplit(graph, lowerBound, upperBound, enforceTrace = True):
    k = len(graph.context.subalignments)
    baseIdx = max(range(k), key = lambda x : upperBound[x] - lowerBound[x])
    #baseIdx = min(range(k), key = lambda x : upperBound[x] - lowerBound[x] if upperBound[x] - lowerBound[x] >= 2 else float('inf'))
    baseLength = upperBound[baseIdx] - lowerBound[baseIdx]
    if baseLength < 2:
        return [list(lowerBound), list(upperBound)]
    
    clusters = initialSplitExpansion(graph, lowerBound, upperBound, baseIdx, baseLength)
    #clusters = initialSplitExpansionSimple(graph, lowerBound, upperBound, baseIdx, baseLength)
    
    cuts = clustersToCuts(graph, lowerBound, upperBound, clusters)   
        
    return cuts

def initialSplitExpansion(graph, lowerBound, upperBound, baseIdx, baseLength):
    k = len(graph.context.subalignments)
    clusters = [[lowerBound[baseIdx] + i] for i in range(baseLength)]
    #idxSets = [set([baseIdx]) for i in range(baseLength)]
    idxSets = {(i, baseIdx) : lowerBound[baseIdx] + i for i in range(baseLength)}
    usedNodes = set()
    weightMap = {}
    
    boundsMap = {}
    for i in range(k):
        boundsMap[0, i] = (lowerBound[i]-1, upperBound[i])
        boundsMap[baseLength-1, i] = (lowerBound[i]-1, upperBound[i])            
    
    heap = []
    for node in range(lowerBound[baseIdx], upperBound[baseIdx]):
        for nbr, value in graph.matrix[node].items():
            i, pos = graph.matSubPosMap[nbr]
    #for i in range(k):
    #    for nbr, value in graph.nodeEdges[node][i]:
            
            if nbr < lowerBound[i]:
                continue
            if nbr >= upperBound[i]:
                continue
                #break
            idx = node - lowerBound[baseIdx]
            if (idx, i) in idxSets:
                continue
            
            heapq.heappush(heap, (-1*value, node, nbr, idx)) 
            weightMap[idx, nbr] = value
    
    while len(heap) > 0:
        value, a, b, idx = heapq.heappop(heap)
        if b in usedNodes:
            continue
        #asub, apos = graph.matSubPosMap[a]
        bsub, bpos = graph.matSubPosMap[b]

        if (idx, bsub) in idxSets:
            continue
        lower, upper = getBounds(boundsMap, baseLength, idx, bsub)
        if not (b > lower and b < upper):
            continue
        
        addBounds(graph, boundsMap, baseLength, idx, b)
        clusters[idx].append(b)
        idxSets[idx, bsub] = b
        usedNodes.add(b)
        
        for nbr, value in graph.matrix[b].items():
            i, pos = graph.matSubPosMap[nbr]
    #for i in range(k):
            if (idx, i) in idxSets:
                continue
            lower, upper = getBounds(boundsMap, baseLength, idx, i)
        #for nbr, value in graph.nodeEdges[b][i]:
            if nbr in usedNodes:
                continue               
             
            if nbr <= lower:
                continue
            if nbr >= upper:
                #break
                continue
            
            #print(weightMap.get((idx, nbr), 0))
            weight = value + weightMap.get((idx, nbr), 0)
            weightMap[idx, nbr] = weight
            heapq.heappush(heap, (-1*weight, b, nbr, idx))    

    return clusters


def getBounds(boundsMap, baseLength, idx, asub):
    a, b = 0, baseLength - 1
    if idx == a:
        return boundsMap[a, asub]
    if idx == b:
        return boundsMap[b, asub]
    
    midpoint = int((a+b)*0.5)
    while (midpoint, asub) in boundsMap:
        if idx == midpoint:
            return boundsMap[midpoint, asub]
        elif idx > midpoint:
            a = midpoint
        elif idx < midpoint:
            b = midpoint
        midpoint = int((a+b)*0.5)
    la, ua = boundsMap[a, asub]
    lb, ub = boundsMap[b, asub]
    return (la, ub)

def addBounds(graph, boundsMap, baseLength, idx, node):
    asub, apos = graph.matSubPosMap[node]
    a, b = 0, baseLength - 1
    
    while True:
        la, ua = boundsMap[a, asub]
        lb, ub = boundsMap[b, asub]
        if idx == a:
            boundsMap[a, asub] = (node, node)
            return
        elif node < ua:
            boundsMap[a, asub] = (la, node)
            
        if idx == b:
            boundsMap[b, asub] = (node, node)
            return
        elif node > lb:
            boundsMap[b, asub] = (node, ub)    
                
        midpoint = int((a+b)*0.5)
        if idx == midpoint:
            boundsMap[midpoint, asub] = (node, node)
            return
        elif (midpoint, asub) not in boundsMap:
            boundsMap[midpoint, asub] = (la, ub)
            
        if idx > midpoint:
            a = midpoint
        elif idx < midpoint:
            b = midpoint 
    

def initialSplitExpansionSimple(graph, lowerBound, upperBound, baseIdx, baseLength):
    k = len(graph.context.subalignments)
    clusters = [[lowerBound[baseIdx] + i] for i in range(baseLength)]
    #idxSets = [set([baseIdx]) for i in range(baseLength)]
    idxSets = {(i, baseIdx) : lowerBound[baseIdx] + i for i in range(baseLength)}
    usedNodes = set()
    weightMap = {}
    
    heap = []
    for node in range(lowerBound[baseIdx], upperBound[baseIdx]):
        for i in range(k):
            for nbr, value in graph.nodeEdges[node][i]:
                if nbr < lowerBound[i]:
                    continue
                if nbr >= upperBound[i]:
                    break
                idx = node - lowerBound[baseIdx]
                heapq.heappush(heap, (-1*value, node, nbr, idx)) 
                weightMap[idx, nbr] = value
    
    while len(heap) > 0:
        value, a, b, idx = heapq.heappop(heap)
        if b in usedNodes:
            continue
        #asub, apos = graph.matSubPosMap[a]
        bsub, bpos = graph.matSubPosMap[b]
        #if bsub in idxSets[idx]:
        #    continue
        if (idx, bsub) in idxSets or idxSets.get((idx-1, bsub), 0) > b or idxSets.get((idx+1, bsub), upperBound[bsub]) < b:
            continue
        
        
        clusters[idx].append(b)
        #idxSets[idx].add(bsub)
        idxSets[idx, bsub] = b
        usedNodes.add(b)
        
        for i in range(k):
            if (idx, i) in idxSets:
                continue
            for nbr, value in graph.nodeEdges[b][i]:
                if nbr in usedNodes:
                    continue                
                if nbr <  idxSets.get((idx-1, i), lowerBound[i]):
                    continue
                if nbr >= idxSets.get((idx+1, i), upperBound[i]):
                    break
                
                #print(weightMap.get((idx, nbr), 0))
                weight = value + weightMap.get((idx, nbr), 0)
                #weight = value
                weightMap[idx, nbr] = weight
                heapq.heappush(heap, (-1*weight, b, nbr, idx))    

    return clusters

def clustersToCuts(graph, lowerBound, upperBound, clusters):
    cuts = [list(lowerBound)]
    #cuts = []
    cut = list(lowerBound)
    for i, cluster in enumerate(clusters):
        if i == 0:
            continue
        #cut = list(cuts[-1])
        
        for a in cluster:
            asub, apos = graph.matSubPosMap[a]
            cut[asub] = max(a, cut[asub])
        cuts.append(cut)
        cut = list(cut)    
    cuts.append(list(upperBound))
    return cuts

def cutsToClusters(graph, cuts):
    clusters = []
    for i in range(len(cuts)-1):
        cluster = []
        for j in range(len(cuts[i])):
            cluster.extend(list(range(cuts[i][j], cuts[i+1][j])))
        clusters.append(cluster)
        #print(cluster)
        
    return clusters