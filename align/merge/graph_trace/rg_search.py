'''
Created on Aug 23, 2020

@author: Vlad
'''

import heapq
from collections import deque 

from configuration import Configs

'''
Region-growing strategy, similar to Kruskal's algorithm. Start with an empty graph, greedily
adding the heaviest edges first. Can be used for clustering or tracing.
'''

def rgSearch(graph):
    Configs.log("Finding graph trace with region-growing search..")
    
    k = len(graph.context.subalignments)
    lowerBound = [graph.subsetMatrixIdx[i] for i in range(k)]
    upperBound = [graph.subsetMatrixIdx[i] + graph.subalignmentLengths[i] for i in range(k)] 
    graph.clusters = rgCluster(graph, lowerBound, upperBound, True)

def rgCluster(graph, lowerBound, upperBound, enforceTrace = True):
    clusters = []
    clusterPointers = {}
    clusterPos = {}
    nodeClusters = {}
    weightMap = []
    absorbed = set()
    cantConnects = set()
    
    for s in range(len(lowerBound)):
        for a in range(lowerBound[s], upperBound[s]):
            clusters.append([a])
            idx = len(clusters)-1
            nodeClusters[a] = idx
            weightMap.append({})
            clusterPos[idx] = {s : a}
            clusterPointers[idx] = {s : (idx-1 if idx > lowerBound[s] else None, idx+1 if idx < upperBound[s]-1 else None)}
    
    
    heap = buildHeap(graph, nodeClusters, weightMap, lowerBound, upperBound)
    Configs.log("Built a heap of size {}..".format(len(heap)))
    crunchHeap(graph, heap, clusters, nodeClusters, clusterPos, clusterPointers, weightMap, cantConnects, absorbed, enforceTrace)

    #c2 = [sorted(c) for c in clusters if len(c) > 0]
    #c2.sort(key= lambda l : graph.matSubPosMap[l[0]])
    #for c in c2:
    #    print(c)
        
    if enforceTrace:
        clusters = orderClusters(graph, clusters, nodeClusters, lowerBound, upperBound)
        #for c in clusters:
        #    print(sorted(c))
    return clusters

def buildHeap(graph, nodeClusters, weightMap, lowerBound, upperBound):
    heap = []
    for s in range(len(lowerBound)):
        for a in range(lowerBound[s], upperBound[s]):
            asub, apos = graph.matSubPosMap[a]
            i = nodeClusters[a]
            for b, value in graph.matrix[a].items():
                bsub, bpos = graph.matSubPosMap[b]
                if b <= a or asub == bsub or b < lowerBound[bsub] or b >= upperBound[bsub]:
                    continue
                j = nodeClusters[b]
                weightMap[j][i] = value
                weightMap[i][j] = value
                heapq.heappush(heap, (-1 * value, a, b))
    
    return heap
    #baseIdx = max(range(k), key = lambda x : upperBound[x] - lowerBound[x])
    #baseLength = upperBound[baseIdx] - lowerBound[baseIdx]

def crunchHeap(graph, heap, clusters, nodeClusters, clusterPos, clusterPointers, weightMap, cantConnects, absorbed, enforceTrace):
    while len(heap) > 0:
        value, a, b = heapq.heappop(heap)
        i, j = nodeClusters[a], nodeClusters[b]
        if i == j or orderPair(i,j) in cantConnects:
            continue
        
        if not checkConnect(graph, i, j, clusters, clusterPos, enforceTrace):
            cantConnects.add(orderPair(i,j))
            continue
        
        absorbed.add(j)
        for e in clusters[j]:
            nodeClusters[e] = i
            clusters[i].append(e)
            asub, apos = graph.matSubPosMap[e]
            clusterPos[i][asub] = e
        clusters[j] = []

        if enforceTrace:        
            for s in clusterPointers[j]:
                prev, nxt = clusterPointers[j][s]
                if prev is not None:
                    clusterPointers[prev][s] = (clusterPointers[prev][s][0], i) 
                if nxt is not None:
                    clusterPointers[nxt][s] = (i, clusterPointers[nxt][s][1])
                clusterPointers[i][s] = (prev, nxt)
    
            updateMergePointers(graph, i, clusterPointers, clusters, clusterPos)
         
        #print("Clusters left: {}".format(len(clusters) - len(absorbed)))   
        for n in weightMap[j]:
            if n in absorbed:
                continue
            weightMap[i][n] = weightMap[i].get(n, 0) + weightMap[j][n]
            weightMap[n][i] = weightMap[i][n]
            heapq.heappush(heap, (-1 * weightMap[i][n], clusters[i][0], clusters[n][0]))

def updateMergePointers(graph, i, clusterPointers, clusters, clusterPos):
    subsets = [graph.matSubPosMap[a][0] for a in clusters[i]]
    
    for s in subsets:
        queue = deque([i])
        visited = set([i])
        
        while len(queue) > 0:
            curNode = queue.popleft()
            
            if clusterPos[curNode].get(s, float('inf')) > clusterPos[i][s] or curNode == i:
                clusterPos[curNode][s] = clusterPos[i][s]
                
                for p in clusterPointers[curNode]:
                    prv, nxt = clusterPointers[curNode][p]
                    if prv not in visited and prv is not None:
                        queue.append(prv)
                        visited.add(prv)
           
            
def checkConnect(graph, i, j, clusters, clusterPos, enforceTrace):
    ci , cj = set([graph.matSubPosMap[a][0] for a in clusters[i]]), set([graph.matSubPosMap[a][0] for a in clusters[j]])
    for s in ci:
        if s in cj:
            return False
    
    if not enforceTrace:
        return True
    
    for s in ci:
        if clusterPos[j].get(s, float('inf')) <= clusterPos[i][s]:
            return False
    for s in cj:
        if clusterPos[i].get(s, float('inf')) <= clusterPos[j][s]:
            return False    
    return True

def orderClusters(graph, clusters, nodeClusters, lowerBound, upperBound):
    orderedClusters = []
    frontier = list(lowerBound)
    while True:
        foundGood = False
        for j in range(len(lowerBound)):
            good = True
            idx = frontier[j]
            if idx >= upperBound[j]:
                continue
            i = nodeClusters[idx]
            for b in clusters[i]:
                bsub, bpos = graph.matSubPosMap[b]
                if b > frontier[bsub]:
                    #print(bsub, b, frontier[bsub])
                    good = False
                    break
                
            if good:
                orderedClusters.append(clusters[i])
                for b in clusters[i]:
                    bsub, bpos = graph.matSubPosMap[b]
                    frontier[bsub] = b + 1
                foundGood = True
                break
        if not foundGood:
            break
    return orderedClusters

def orderPair(a, b):
    return (min(a, b), max(a, b))