'''
Created on Jun 15, 2020

@author: Vlad
'''

import heapq
import time 
import random

from magus_configuration import Configs
    
'''
Fiduccia-Mattheyses implementation for clustering and/or tracing. 
Currently not recommended, slower and less accurate than minclusters or mwtgreedy.
'''

def fmAlgorithm(graph):
    Configs.log("Finding graph trace with FM Algorithm..")
    
    k = len(graph.context.subalignments)
    lowerBound = [graph.subsetMatrixIdx[i] for i in range(k)]
    upperBound = [graph.subsetMatrixIdx[i] + graph.subalignmentLengths[i] for i in range(k)]  
    
    if graph.clusters is None or len(graph.clusters) == 0:
        graph.buildNodeEdgeDataStructure()
    else:
        graph.buildNodeEdgeDataStructureFromClusters()
    clusters, totalCost, cuts = fmPartition(graph, lowerBound, upperBound)
    
    graph.clusters = clusters
                    
def fmPartitionWithCuts(graph, lowerBound, upperBound, cuts):
    allCuts = [lowerBound] + cuts + [upperBound]
    clusters, cost, cuts = [], 0, []
    for i in range(len(allCuts)-1):
        newClusters, newCost, newCuts = fmPartition(graph, allCuts[i], allCuts[i+1], False)
        clusters = clusters + newClusters
        cost = cost + newCost
        cuts = cuts + newCuts    
    return clusters, cost, cuts
    
def fmPartition(graph, lowerBound, upperBound, iterate = True):
    k = len(graph.context.subalignments)
    
    finished = True
    cluster = []
    for i in range(k):
        if upperBound[i] - lowerBound[i] > 1:
            finished = False
            break
        elif upperBound[i] - lowerBound[i] == 1:
            cluster.append(lowerBound[i])
    if finished:
        return [cluster], 0, []   
    
    startingCut = [int((lowerBound[i] + upperBound[i])*0.5) for i in range(k)]    
    bestCut, bestCutCost = fmFindBestCut(graph, lowerBound, upperBound, startingCut, None)    
        
    lowerClusters, lowerCost, lowerCuts = fmPartition(graph, lowerBound, bestCut)
    upperClusters, upperCost, upperCuts = fmPartition(graph, bestCut, upperBound)
    clusters = lowerClusters + upperClusters
    totalCost = bestCutCost + lowerCost + upperCost
    totalCuts = lowerCuts + [bestCut] + upperCuts

    return clusters, totalCost, totalCuts

def fmFindBestCut(graph, lowerBound, upperBound, startingCut, widthSumLimit):
    #Configs.log("Finding FM partition..")
    #Configs.log("    Lower bound: {}".format(graph.cutString(lowerBound)))
    #Configs.log("    Upper bound: {}".format(graph.cutString(upperBound)))
    #Configs.log("    Startng cut: {}".format(graph.cutString(startingCut)))
    k = len(graph.context.subalignments)
    bestCut = startingCut
    bestCutGains, bestCutCost = populateGains(graph, lowerBound, upperBound, bestCut)
    
    oldBestCut = None
    while bestCut != oldBestCut:

        if bestCutCost == 0:
            break
        
        oldBestCut = bestCut
        cutCost = bestCutCost
        cut = list(bestCut)
        gains = dict(bestCutGains)
        heap = []
        heapGains = {}
        heapGainsVersions = {}
        locked = set()
        lowerUpdateBound = [cut[i] for i in range(k)]
        upperUpdateBound = [cut[i] for i in range(k)]
        
        while True:
            newLowerBound, newUpperBound = findNewBounds(graph, lowerBound, upperBound, cut, widthSumLimit)            
            updateList = getHeapGainUpdateList(graph, newLowerBound, newUpperBound, lowerUpdateBound, upperUpdateBound)   
            updateHeapGainList(graph, cut, updateList, gains, heapGains, heapGainsVersions, heap, locked) 
            
            found = False
            reinsert = []
            while len(heap) > 0:
                gain, node, gainVersion = heapq.heappop(heap)
                if node in locked or gainVersion != heapGainsVersions[node]:
                    continue
                
                asub, apos = graph.matSubPosMap[node]
                if node < newLowerBound[asub] or node >= newUpperBound[asub]:
                    reinsert.append((gain, node, gainVersion))
                    continue

                found = True
                break
            
            if not found:
                break
            for item in reinsert:
                heapq.heappush(heap, item)
            
            locked.add(node)
            gain = gain * -1
            #asub, apos = graph.matSubPosMap[node]
            movedNodes = []
            if node >= cut[asub]:
                movedNodes = [j for j in range(cut[asub], node+1)]
                cut[asub] = node+1
            else:
                movedNodes = [j for j in range(node, cut[asub])]
                cut[asub] = node
            
            updateGains(graph, lowerBound, upperBound, cut, gains, movedNodes, lowerUpdateBound, upperUpdateBound) 
            
            cutCost = cutCost - gain
            if cutCost < bestCutCost:
                bestCut = list(cut)
                bestCutGains = dict(gains)
                bestCutCost = cutCost
                
    Configs.log("Found FM partition {}".format(graph.cutString(bestCut)))
    Configs.log("    Partition cost: {}".format(bestCutCost))  
    return bestCut, bestCutCost

def populateGains(graph, lowerBound, upperBound, cut):
    k = len(graph.context.subalignments)
    gains = {}
    cutCost = 0
    
    for j in range(k):
        for node in range(lowerBound[j], upperBound[j]):
            gain = 0
            for i in range(k):
                for nbr, value in graph.nodeEdges[node][i]:
                    if nbr < lowerBound[i]:
                        continue
                    if nbr >= upperBound[i]:
                        break
                    
                    if (nbr < cut[i] and node < cut[j]) or (nbr >= cut[i] and node >= cut[j]):
                        gain = gain - value
                    else:
                        gain = gain + value
                        cutCost = cutCost + value
            gains[node] = gain 
    return gains, int(cutCost/2)               

def findNewBounds(graph, lowerBound, upperBound, cut, widthSumLimit):
    k = len(graph.context.subalignments)
    lowerSize = sum([cut[i] - lowerBound[i] for i in range(k)])
    upperSize = sum([upperBound[i] - cut[i] for i in range(k)])
    portionWidth = getPortionWidth(lowerBound, upperBound)
    
    limit = min(k, lowerSize + upperSize - 1)
    
    lowerMargin = int((lowerSize - upperSize + limit) * 0.5)
    upperMargin = int((upperSize - lowerSize + limit) * 0.5)
    
    #newLowerBound = [max(cut[i]-lowerMargin, lowerBound[i]) for i in range(k)]
    #newUpperBound = [min(cut[i]+upperMargin, upperBound[i]) for i in range(k)]
    newLowerBound = [max(cut[i]-lowerMargin, lowerBound[i], upperBound[i] - portionWidth + 1) for i in range(k)]
    newUpperBound = [min(cut[i]+upperMargin, upperBound[i], lowerBound[i] + portionWidth - 1) for i in range(k)]
    
    if widthSumLimit is not None:
        l1, l2, u1, u2 = None, None, None, None
        for i in range(k):
            l, u = cut[i]-lowerBound[i], upperBound[i]-cut[i]
            if l1 is None or l > cut[l1]-lowerBound[l1]:
                l1, l2 = i, l1
            elif l2 is None or l > cut[l2]-lowerBound[l2]:
                l2 = i
            if u1 is None or u > upperBound[u1]-cut[u1]:
                u1, u2 = i, u1
            elif u2 is None or u > upperBound[u2]-cut[u2]:
                u2 = i
        
        for i in range(k):
            if i == l1:
                newLowerBound[i] = max(newLowerBound[i], cut[l2] - lowerBound[l2] + upperBound[i] - widthSumLimit)
            else:
                newLowerBound[i] = max(newLowerBound[i], cut[l1] - lowerBound[l1] + upperBound[i] - widthSumLimit)
            if i == u1:
                newUpperBound[i] = min(newUpperBound[i], widthSumLimit + lowerBound[i] + cut[u2] - upperBound[u2])
            else:
                newUpperBound[i] = min(newUpperBound[i], widthSumLimit + lowerBound[i] + cut[u1] - upperBound[u1])
    
    return newLowerBound, newUpperBound

def updateGains(graph, lowerBound, upperBound, cut, gains, movedNodes, lowerUpdateBound, upperUpdateBound):
    k = len(graph.context.subalignments)
        
    for node in movedNodes:
        gains[node] = -1 * gains[node]
        asub, apos = graph.matSubPosMap[node]
        lowerUpdateBound[asub] = cut[asub]
        upperUpdateBound[asub] = cut[asub]
        
        for i in range(k):
            for nbr, value in graph.nodeEdges[node][i]:
                if nbr < lowerBound[i]:
                    continue
                if nbr >= upperBound[i]:
                    break
                
                if (nbr < cut[i] and node < cut[asub]) or (nbr >= cut[i] and node >= cut[asub]):
                    gains[nbr] = gains[nbr] - 2 * value
                else:
                    gains[nbr] = gains[nbr] + 2 * value
                
                if nbr < cut[i]:
                    lowerUpdateBound[i] = max(lowerUpdateBound[i], nbr+1)
                else:
                    upperUpdateBound[i] = min(upperUpdateBound[i], nbr)

def getHeapGainUpdateList(graph, newLowerBound, newUpperBound, lowerUpdateBound, upperUpdateBound):
    k = len(graph.context.subalignments)
    updateList = [[] for i in range(k)]      
    
    for i in range(k):
        for j in reversed(range(newLowerBound[i], lowerUpdateBound[i])):
            updateList[i].append(j)
        
        for j in range(upperUpdateBound[i], newUpperBound[i]):  
            updateList[i].append(j)
    
    return updateList

def updateHeapGainList(graph, cut, updateList, gains, heapGains, heapGainsVersions, heap, locked):
    k = len(graph.context.subalignments)      
    
    for i in range(k):
        for j in updateList[i]:
            heapGainsVersions[j] = heapGainsVersions.get(j, 0) + 1
            
            if j > cut[i]:
                heapGains[j] = heapGains[j-1] + gains[j]
            elif j == cut[i] or j == cut[i]-1:
                heapGains[j] = gains[j]
            else:
                heapGains[j] = heapGains[j+1] + gains[j]
            
            if j not in locked:
                heapq.heappush(heap, (-heapGains[j], j, heapGainsVersions[j]))

def getPortionWidth(lowerBound, upperBound):
    return max([u-l for u,l in zip(upperBound, lowerBound)])

def computeCutCost(graph, lowerBound, upperBound, cut):
    k = len(graph.context.subalignments)
    cutCost = 0
    
    for j in range(k):
        for node in range(lowerBound[j], upperBound[j]):
            for i in range(j, k):
                for nbr, value in graph.nodeEdges[node][i]:
                    if nbr < lowerBound[i]:
                        continue
                    if nbr >= upperBound[i]:
                        break
                    
                    if (nbr >= cut[i] and node < cut[j]) or (nbr < cut[i] and node >= cut[j]):
                        cutCost = cutCost + value
    return cutCost 