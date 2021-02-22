'''
Created on Apr 18, 2020

@author: Vlad
'''

import heapq

from configuration import Configs

'''
Resolve clusters into a trace by breaking conflicting clusters apart.
We use A* to search for the path of cluster breaks with the smallest number of clusters broken.
'''

#todo refactoring
def minClustersSearch(graph):
    Configs.log("Finding graph trace with minimum clusters heuristic search..")
    
    subsetClusters = {}
    clusterPositions = {}
    queueIdxs = {}
    clusterBreaks = {}
    maxFrontier = {}
    visitedStates = set()
    maximalCut = {}
    stateCounter = 0
    aggression = 1.0
    lastFrontierState = None
    greedy = False
    totalPairs = 0

    
    for a,cluster in enumerate(graph.clusters):
        for b in cluster:
            bsub, bpos = graph.matSubPosMap[b] 
            subsetClusters[bsub] = subsetClusters.get(bsub, []) + [(a, bpos)]  
            clusterPositions[a] = {}
            totalPairs = totalPairs + len(cluster)*(len(cluster)-1)/2
                
    
    for asub in subsetClusters:
        subsetClusters[asub].sort(key = lambda c: c[1])
        queueIdxs[asub] = 0

        maxFrontier[asub] = -1
        maximalCut[asub] = -1
        
        for i in range(len(subsetClusters[asub])):            
            a = subsetClusters[asub][i][0]
            clusterPositions[a][asub] = i
    
    heap = []
    
    startState = (0, 0, len(graph.clusters), totalPairs, stateCounter, queueIdxs, clusterBreaks, maximalCut, [], True)
    startState = developState(startState, graph, aggression, greedy, 0, subsetClusters, clusterPositions)
    heapq.heappush(heap, startState)
    
    while len(heap) > 0:
        heapCleared = False
        if len(heap) > Configs.searchHeapLimit:
            Configs.log("Heap limit {} reached.. Truncating heap to last frontier".format(Configs.searchHeapLimit)) 
            if aggression == 1:
                aggression = 1.2
                Configs.log("Increasing aggression to {}..".format(aggression))
            elif aggression < 8:
                aggression = int(aggression) * 2
                Configs.log("Increasing aggression to {}..".format(aggression))
            else:
                Configs.log("Setting search strategy to fully greedy..")
                greedy = True
                aggression = 1    
                
            heap.clear()
            visitedStates = set()
            lastFrontierState = developState(lastFrontierState, graph, aggression, greedy, 0, subsetClusters, clusterPositions)
            heapq.heappush(heap, lastFrontierState)
            heapCleared = True

        state = heapq.heappop(heap)
        heuristic, numOrdered, numLeft, pairsLeft, counter, queueIdxs, clusterBreaks, maximalCut, newClusterBreaks, safeFrontier = state
   
        if len(newClusterBreaks) == 0:
            break
        else:
            stateKey = tuple([queueIdxs[sub] for sub in subsetClusters])
            if stateKey in visitedStates:
                continue
            else:
                visitedStates.add(stateKey)
            

            newSubFrontier = True
            for asub in queueIdxs:
                if queueIdxs[asub] <= maxFrontier[asub]:
                    newSubFrontier = False
                    break
                
            if newSubFrontier:
                maxFrontier = queueIdxs
                Configs.log("Reached new search frontier")
                Configs.log(maxFrontier)
                lastFrontierState = state
                greedy = False
            
            if safeFrontier and not heapCleared:
                Configs.log("Safe frontier reached.. dumping {} from heap and resetting aggression..".format(len(heap)))
                lastFrontierState = state
                heap.clear()
                visitedStates = set()
                aggression = 1.0
                greedy = False
                
            
            nextStates = []
            for a, goodSide, badSide, crossedClusters in newClusterBreaks:
                
                g,b = len(goodSide), len(badSide)
                pairsDiff = g*(g-1)/2 + b*(b-1)/2 - (g+b)*(g+b-1)/2
                
                stateCounter = stateCounter + 1
                queueIdxsCopy = dict(queueIdxs)
                clusterBreaksCopy = dict(clusterBreaks)
                maximalCutCopy = dict(maximalCut)
                for b in goodSide:
                    bsub, bpos = graph.matSubPosMap[b]
                    clusterBreaksCopy[a, bsub] = goodSide
                    maximalCutCopy[bsub] = max(maximalCutCopy[bsub], clusterPositions[a][bsub])
                    
                for b in badSide:
                    bsub, bpos = graph.matSubPosMap[b]
                    clusterBreaksCopy[a, bsub] = badSide
                    maximalCutCopy[bsub] = max(maximalCutCopy[bsub], clusterPositions[a][bsub])
                
                nextState = (0, numOrdered, numLeft + 1, pairsLeft + pairsDiff, stateCounter, queueIdxsCopy, clusterBreaksCopy, maximalCutCopy, [], False)
                nextState = developState(nextState, graph, aggression, greedy, len(crossedClusters), subsetClusters, clusterPositions)

                nextStates.append(nextState)
                
            if greedy:
                nextState = min(nextStates, key=lambda x : x[0])
                heapq.heappush(heap, nextState)
            else:
                for nextState in nextStates:
                    heapq.heappush(heap, nextState)
                
    queueIdxs = {}
    for asub in subsetClusters:
        queueIdxs[asub] = 0
        
    orderedClusters = []
    foundGood = True
    while foundGood:
        foundGood = False
        for asub in queueIdxs:
            good = True
            idx = queueIdxs[asub]
            if idx == len(subsetClusters[asub]):
                continue
            a, pos = subsetClusters[asub][idx]
            
            if (a,asub) in clusterBreaks:
                cluster = clusterBreaks[a,asub]
            else:
                cluster = graph.clusters[a]
            
            for b in cluster:
                bsub, bpos = graph.matSubPosMap[b]
                if clusterPositions[a][bsub] != queueIdxs[bsub]:
                    good = False
                    break
                
            if good:
                orderedClusters.append(cluster)
                for b in cluster:
                    bsub, bpos = graph.matSubPosMap[b]
                    queueIdxs[bsub] = clusterPositions[a][bsub] + 1
                foundGood = True
                break
                    
    graph.clusters = orderedClusters

def developState(state, graph, aggression, greedy, crossed, subsetClusters, clusterPositions):
    heuristic, numOrdered, numLeft, pairsLeft, counter, queueIdxs, clusterBreaks, maximalCut, newClusterBreaks, safeFrontier  = state
    
    foundGood = True
    while foundGood:
        foundGood = False
        newClusterBreaks = []
        visited = set()
        safeFrontier = True
        
        for asub in queueIdxs:
            idx = queueIdxs[asub]
            if idx <= maximalCut[asub]:
                safeFrontier = False

            if idx == len(subsetClusters[asub]):
                continue
            a, pos = subsetClusters[asub][idx]
            if (a,asub) in visited:
                continue
            
            if (a,asub) in clusterBreaks:
                cluster = clusterBreaks[a,asub]
            else:
                cluster = graph.clusters[a]
            
            goodSide, badSide, crossedClusters = [], [], set()
            for b in cluster:
                bsub, bpos = graph.matSubPosMap[b]
                visited.add((a, bsub))
                bidx = clusterPositions[a][bsub]
                diff = bidx - queueIdxs[bsub]
                if diff == 0:
                    goodSide.append(b)
                else:
                    badSide.append(b)
                                   
            if len(badSide) == 0:
                for b in cluster:
                    bsub, bpos = graph.matSubPosMap[b]
                    queueIdxs[bsub] = clusterPositions[a][bsub] + 1
                numOrdered = numOrdered + 1
                numLeft = numLeft - 1
                foundGood = True
                break
            else:
                newClusterBreaks.append((a, goodSide, badSide, crossedClusters))

    
    if greedy:
        for a, goodSide, badSide, crossedClusters in newClusterBreaks:
            goodSub = set()
            for b in goodSide:
                bsub, bpos = graph.matSubPosMap[b] 
                goodSub.add(bsub)
                
            for b in badSide:    
                bsub, bpos = graph.matSubPosMap[b]   
                for i in range(queueIdxs[bsub], clusterPositions[a][bsub]):
                    c, posc = subsetClusters[bsub][i]
                    if (c,bsub) in clusterBreaks:
                        otherCluster = clusterBreaks[c,bsub]
                    else:
                        otherCluster = graph.clusters[c]
                    
                    for csite in otherCluster:
                        csub, cpos = graph.matSubPosMap[csite]
                        if csub in goodSub and clusterPositions[c][csub] > clusterPositions[a][csub]:
                            crossedClusters.add(c)
                            break  
                          
                
    if safeFrontier or len(newClusterBreaks) == 0:
        heuristic = (numLeft + numOrdered, -numOrdered, -crossed, -pairsLeft)
    else:
        heuristic = (aggression * numLeft + numOrdered, -numOrdered, -crossed, -pairsLeft)
    state = (heuristic, numOrdered, numLeft, pairsLeft, counter, queueIdxs, clusterBreaks, maximalCut, newClusterBreaks, safeFrontier)
    return state