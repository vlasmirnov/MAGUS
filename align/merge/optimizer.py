'''
Created on Aug 23, 2020

@author: Vlad
'''

import heapq
import time

from configuration import Configs

'''
Optimizer may be used to post-process a trace to improve the MWT score by shuffling nodes between clusters.
Disabled by default - tends to be very time-consuming with negligible improvements to accuracy.
However, may be helpful when using inaccurate clustering and/or tracing algorithms.
The quality of the resulting trace is usually on par with using MCL/minclusters.
'''

def optimizeTrace(graph):
    time1 = time.time() 
    
    if Configs.graphTraceOptimize:
        Configs.log("Optimization pass..")
        graph.addSingletonClusters()
        graph.clusters = optimizeClusters(graph, graph.clusters)
        Configs.log("Optimized the trace to {} clusters with a total cost of {}".format(len(graph.clusters), graph.computeClusteringCost(graph.clusters)))
    else:
        Configs.log("Skipping optimization pass..")
    time2 = time.time()
    Configs.log("Finished optimization in {} sec..".format(time2-time1))
    

def optimizeClusters(graph, clusters):
    bestClusters, bestCost = clusters, graph.computeClusteringCost(clusters)
    Configs.log("Starting optimization from initial cost of {}..".format(bestCost)) 
    context = SearchContext(clusters)
    context.initialize(graph)
    
    passNum = 1
    while True:
        Configs.log("Starting optimization pass {}..".format(passNum))
        newClusters, gain = optimizationPass(graph, bestClusters, context)
        if gain > 0:
            bestClusters = newClusters
            bestCost = bestCost - gain
            Configs.log("New clustering with a cost of {} over {} clusters..".format(bestCost, len(bestClusters)))
            #Configs.log("Verifying cost of {}..".format(graph.computeClusteringCost(bestClusters)))
        else:
            break
        passNum = passNum + 1
    #Configs.log("Final optimized cost of {} over {} clusters..".format(graph.computeClusteringCost(bestClusters), len(bestClusters)))
    Configs.log("Final optimized cost of {} over {} clusters..".format(bestCost, len(bestClusters)))
    return bestClusters

def optimizationPass(graph, clusters, context):
    context.initializeHeap(graph)
    #cost = graph.computeClusteringCost(clusters)
    #Configs.log("Initial cost: {}..".format(cost))
    
    bestGain, currentGain, bestClusters = 0, 0, clusters
    
    #move = 0
    while True:
        nextMove = getNextClusterMove(graph, context)
        if nextMove is None:
            if currentGain == bestGain:
                bestClusters = context.getOrderedClusters()
            break
            
        gain, element, src, dest, updateList = nextMove        
        if gain < 0 and currentGain == bestGain:
            bestClusters = context.getOrderedClusters()
            break
        
        currentGain = currentGain + gain
        if currentGain > bestGain:
            bestGain = currentGain
            context.locked = set()
        else:
            context.locked.add(element)
            
        context.moveElements(graph, element, dest, updateList)
        context.updateMoves(graph, element, src, dest, updateList) 
        
        #move = move + 1
        #if move % 10000 == 0:       
        #Configs.log("Current / best cost: {} / {}..".format(cost - currentGain, cost - bestGain))
    return bestClusters, bestGain

def getNextClusterMove(graph, context):
    while len(context.heap) > 0:
        gain, element, dest = heapq.heappop(context.heap)
        gain = gain * -1
        src = context.elementClusters[element]
        context.elementMoves[element, dest] = None
        
        if element in context.locked or dest in context.deletedClusters:
            continue
        
        simpleGain = context.getGainSimple(element, dest)
        if simpleGain <= 0:
            continue
        
        updateList = context.getElementUpdateList(graph, element, dest)
        updatedGain = context.getGainCorrectionFast(simpleGain, updateList)
        #Configs.log("Gain is {}, updated gain is {}, heap size {}..".format(gain, updatedGain, len(context.heap)))
        
        if updatedGain < gain:
            if updatedGain > 0:
                context.elementMoves[element, dest] = updatedGain
                heapq.heappush(context.heap, (-1*updatedGain, element, dest))
            continue
        #Configs.log("Completed the move..")
        
        return updatedGain, element, src, dest, updateList
        
    return None      


class SearchContext:
    
    def __init__(self, clusters):
        self.clusters = [set(c) for c in clusters]
        self.clusterLL = {}
        self.clusterLLHead = 0
        self.clusterOrders = {}
        self.elementClusters = {}
        self.clusterSubs = {}
        self.deletedClusters = set()
        
        self.weights = {}
        self.gainStructure = []
        self.elementMoves = {}
        self.heap = []
        self.locked = set()
        
        self.mode = "positive_moves"
        #self.mode = "adjacent_moves"
            
    def initialize(self, graph):
        Configs.log("Initializing search context data structures..")
        for i in range(len(self.clusters)):
            self.clusterOrders[i] = [i]
            self.clusterLL[i] = (i-1 if i > 0 else None, i+1 if i < len(self.clusters)-1 else None)
            for a in self.clusters[i]:
                asub, apos = graph.matSubPosMap[a]
                self.elementClusters[a] = i
                self.clusterSubs[i, asub] = a
                nbrs = self.getNeighborList(graph, a)
                self.updateNeighborWeights(None, i, nbrs)
        #Configs.log("Initialized..")
        
    def initializeHeap(self, graph):
        Configs.log("Reinitializing heap and all that stuff..")
        self.gainStructure = []
        self.elementMoves = {}
        self.heap = []
        self.locked = set()
                
        print("Working with {} clusters..".format(len(self.clusters) - len(self.deletedClusters)))
        
        k = len(graph.context.subalignments)
        self.gainStructure = [[0 for j in range(graph.subalignmentLengths[i])] for i in range(k)]        
        for i in range(k):
            for j in range(graph.subalignmentLengths[i]):
                node = graph.subsetMatrixIdx[i] + j
                weight = self.weights.get((node, self.elementClusters[node]), 0)
                self.gainStructure[i][j] = weight if j == 0 else weight + self.gainStructure[i][j-1]
        
        #if self.mode == "positive_moves":
        #    candidates = self.getPositiveMoves(graph)
        #elif self.mode == "adjacent_moves":
        #    candidates = self.getAdjacentMoves(graph)
                            
        #Configs.log("Considering {} candidate moves..".format(len(candidates)))
        #Configs.log("Choosing {} out of {} candidates..".format(limit, len(candidates)))
        #candidates = heapq.nlargest(limit, candidates)
        #for gain, i, nbr in candidates:
        #    self.pullNeighborMoves(graph, i, [(nbr, gain)])
            #gain = self.getGainSimple(nbr, i)
            #self.elementMoves[nbr, i] = gain
            #heapq.heappush(self.heap, (-1*gain, nbr, i))
        
        self.getPositiveMoves(graph)
        Configs.log("Starting with {} candidate moves..".format(len(self.heap)))
        #import sys
        #sys.exit()
    
    def getPositiveMoves(self, graph):
        k = len(graph.context.subalignments)
        candidates = []
        used = set()        
        
        clusterSubMap = [[-1 for c in self.clusters] for i in range(k)] 
        
        i = self.clusterLLHead
        while i is not None:
            prev, nxt = self.clusterLL[i]
        #for i in range(len(self.clusters)):
            #if i > 0:
            if prev is not None:
                for j in range(k):
                    clusterSubMap[j][i] = clusterSubMap[j][prev]
                    #clusterSubMap[j][i] = clusterSubMap[j][i]
            for node in self.clusters[i]:
                asub, apos = graph.matSubPosMap[node]
                clusterSubMap[asub][i] = apos
            i = nxt
        
        i = self.clusterLLHead
        while i is not None:
            prev, nxt = self.clusterLL[i]                     
        #for i in range(len(self.clusters)):
            
            for node in self.clusters[i]:
                asub, apos = graph.matSubPosMap[node]
                for nbr, value in graph.matrix[node].items():
                    bsub, bpos = graph.matSubPosMap[nbr]
                    j = self.elementClusters[nbr]
                    #if asub == bsub or j <= i:
                    if asub == bsub or self.clusterOrders[j] <= self.clusterOrders[i]:
                        continue
                    
                    if (i, nbr) not in used:
                        used.add((i, nbr))
                        gain = self.weights.get((nbr, i), 0) - self.gainStructure[bsub][bpos]
                        bound = clusterSubMap[bsub][i] - 1 if self.clusterSubs.get((i, bsub)) is not None else clusterSubMap[bsub][i]
                        if bound >= 0:
                            gain = gain + self.gainStructure[bsub][bound]
                        #if gain >= 0:
                        if gain > 0 or (gain == 0 and len(self.clusters[i]) >= len(self.clusters[j])):
                            #candidates.append((gain, i, nbr))
                            self.elementMoves[nbr, i] = gain
                            heapq.heappush(self.heap, (-1*gain, nbr, i))
                        
                    if (j, node) not in used:
                        used.add((j, node))
                        gain = self.weights.get((node, j), 0) - self.gainStructure[asub][clusterSubMap[asub][j]]
                        if apos > 0:
                            gain = gain + self.gainStructure[asub][apos-1]
                        #if gain >= 0:
                        if gain > 0 or (gain == 0 and len(self.clusters[j]) >= len(self.clusters[i])):
                            #candidates.append((gain, j, node))
                            self.elementMoves[node, j] = gain
                            heapq.heappush(self.heap, (-1*gain, node, j))             
            i = nxt
                   
        return candidates
    
    def getAdjacentMoves(self, graph):
        candidates = []
        for i in range(len(self.clusters)-1):
            for a in self.clusters[i]:
                candidates.append((0, i+1, a))
                asub, apos = graph.matSubPosMap[a]
                if a < graph.subsetMatrixIdx[asub] + graph.subalignmentLengths[asub] - 1:
                    for j in range(i, self.elementClusters[a+1]):
                        candidates.append((0, j, a+1))
                
            for a in self.clusters[i+1]:
                candidates.append((0, i, a)) 
                asub, apos = graph.matSubPosMap[a]
                if a > graph.subsetMatrixIdx[asub]:
                    for j in range(self.elementClusters[a-1]+1, i+2):
                        candidates.append((0, j, a-1))       
        return candidates
        
    def getGain(self, element, dest, updateList):
        gain = 0
        for item in updateList:
            gain = gain - self.weights.get((item, self.elementClusters[item]), 0)
        gain = gain - self.weights.get((element, self.elementClusters[element]), 0)
        gain = gain + self.weights.get((element, dest), 0)
        return gain
    
    def getGainSimple(self, element, dest):
        return self.weights.get((element, dest), 0) - self.weights.get((element, self.elementClusters[element]), 0)
    
    def getGainCorrection(self, updateList):
        gain = 0
        for item in updateList:
            gain = gain - self.weights.get((item, self.elementClusters[item]), 0)
        return gain
    
    def getGainCorrectionFast(self, gain, updateList):
        for item in updateList:
            gain = gain - self.weights.get((item, self.elementClusters[item]), 0)
            if gain < 0:
                return gain
        return gain
    
    def getElementUpdateList(self, graph, element, dest):
        asub, apos = graph.matSubPosMap[element]
        updateList = []
        
        curCluster = self.elementClusters[element]
        while curCluster != dest:
            if self.compareClusterOrder(curCluster, dest):
                prev, curCluster = self.clusterLL[curCluster]
            else:
                curCluster, nxt = self.clusterLL[curCluster]
            if self.clusterSubs.get((curCluster, asub)) is not None:
                updateList.append(self.clusterSubs[curCluster, asub])
        return updateList

    def getNeighborList(self, graph, element):
        asub, apos = graph.matSubPosMap[element]
        nbrs = []
        
        for nbr, value in graph.matrix[element].items():
            bsub, bpos = graph.matSubPosMap[nbr]
            if asub != bsub:
                nbrs.append((nbr, value))
        return nbrs
        
    def updateNeighborWeights(self, src, dest, nbrs):
        for nbr, value in nbrs:
            if src is not None:
                self.weights[nbr, src] = self.weights.get((nbr, src), 0) - value
            if dest is not None:
                self.weights[nbr, dest] = self.weights.get((nbr, dest), 0) + value
    
    def pullNeighborMoves(self, graph, dest, nbrs):
        if dest in self.deletedClusters:
            return
        for nbr, value in nbrs:
            if nbr in self.locked or self.elementClusters[nbr] == dest:
                continue               
            
            gain = self.getGainSimple(nbr, dest)
            if self.elementMoves.get((nbr, dest)) is None or gain > self.elementMoves[nbr, dest]:
                #self.elementMoves[nbr, dest] = gain
                #heapq.heappush(self.heap, (-1*gain, nbr, dest))    
                updateList = self.getElementUpdateList(graph, nbr, dest)
                #gain = gain + self.getGainCorrection(updateList)
                gain = self.getGainCorrectionFast(gain, updateList)
                if gain >= 0 and (self.elementMoves.get((nbr, dest)) is None or gain > self.elementMoves[nbr, dest]):
                #if(gain > 0 or (gain == 0 and len(self.clusters[dest]) >= len(self.clusters[self.elementClusters[nbr]]))) \
                #and (self.elementMoves.get((nbr, dest)) is None or gain > self.elementMoves[nbr, dest]):
                    self.elementMoves[nbr, dest] = gain
                    heapq.heappush(self.heap, (-1*gain, nbr, dest))
                #else:
                #    print(nbr, dest, value, gain)

    def moveElement(self, graph, element, src, dest):
        asub, apos = graph.matSubPosMap[element]
        self.clusters[src].remove(element)
        self.clusters[dest].add(element)
        self.elementClusters[element] = dest
        if self.clusterSubs[src, asub] == element:
            self.clusterSubs[src, asub] = None
        self.clusterSubs[dest, asub] = element
        if len(self.clusters[src]) == 0:# and delete:
            self.deleteCluster(src)
        nbrs = self.getNeighborList(graph, element)
        self.updateNeighborWeights(src, dest, nbrs)    
        
    def moveElements(self, graph, element, dest, updateList):
        src = self.elementClusters[element]
        curNode = dest
        self.moveElement(graph, element, src, dest)
        for node in updateList:
            prev, nxt = self.clusterLL[curNode]
            if self.clusterOrders[dest] > self.clusterOrders[src]:
                idx = self.insertCluster(curNode, nxt)
            else:
                idx = self.insertCluster(prev, curNode)
            #self.moveElement(graph, node, self.elementClusters[node], idx, self.elementClusters[node] != dest)
            self.moveElement(graph, node, self.elementClusters[node], idx)
            curNode = idx
        
    def updateMoves(self, graph, element, src, dest, updateList):
        if self.mode == "positive_moves":
            nbrs = self.getNeighborList(graph, element)
            nbrs = [item for item in nbrs if self.getGainSimple(item[0], dest) >= 0]
            self.pullNeighborMoves(graph, dest, nbrs)
            #self.pullNeighborMoves(graph, src, [(element, 0)])
            for node in updateList:
                nbrs = self.getNeighborList(graph, node)
                nbrs = [item for item in nbrs if self.getGainSimple(item[0], self.elementClusters[node]) >= 0]
                self.pullNeighborMoves(graph, self.elementClusters[node], nbrs)
        elif self.mode == "adjacent_moves":
            for a in updateList:
                nbrs = []
                asub, apos = graph.matSubPosMap[a]
                i = self.elementClusters[a]
                prv, nxt = self.clusterLL[i]
                if prv is not None:
                    for b in self.clusters[prv]:
                        nbrs.append((b, 0))
                if nxt is not None:
                    for b in self.clusters[nxt]:
                        nbrs.append((b, 0))
                if a < graph.subsetMatrixIdx[asub] + graph.subalignmentLengths[asub] - 1:
                    nbrs.append((a+1, 0))
                if a > graph.subsetMatrixIdx[asub]:
                    nbrs.append((a-1, 0))
                self.pullNeighborMoves(graph, i, nbrs)
    
    def insertCluster(self, prev, nxt):
        self.clusters.append(set())
        idx = len(self.clusters)-1
        self.clusterOrders[idx] = self.getMiddleOrder(prev, nxt)        
        self.clusterLL[idx] = (prev, nxt)
        if prev is not None:
            self.clusterLL[prev] = (self.clusterLL[prev][0], idx) 
        if nxt is not None:
            self.clusterLL[nxt] = (idx, self.clusterLL[nxt][1])
            
        if prev is None:
            self.clusterLLHead = idx    
        return idx        
    
    def deleteCluster(self, cluster):
        self.deletedClusters.add(cluster)
        prev, nxt = self.clusterLL[cluster]
        if prev is not None:
            self.clusterLL[prev] = (self.clusterLL[prev][0], nxt) 
        if nxt is not None:
            self.clusterLL[nxt] = (prev, self.clusterLL[nxt][1])
        
        if prev is None:
            self.clusterLLHead = nxt
    
    def getMiddleOrder(self, a, b):
        #return (self.clusterOrders[a] + self.clusterOrders[b]) * 0.5
        oa, ob = self.clusterOrders.get(a, []), self.clusterOrders.get(b, [])
        idx = len(oa)-1
        order = list(oa)
        if len(ob) == len(oa):
            order.append(0)
        elif len(ob) < len(oa):            
            order[-1] = order[-1] + 1
        elif len(ob) > len(oa):
            order.append(ob[idx+1]-1)
        
        #print(oa, ob, order)    
        return order
                
    
    def compareClusterOrder(self, a, b):
        oa, ob = self.clusterOrders[a], self.clusterOrders[b]
        idx = 0
        while True:            
            x, y = oa[idx] if idx < len(oa) else -float('inf'), ob[idx] if idx < len(ob) else -float('inf')
            if x < y:
                return True
            if y < x:
                return False
            idx = idx + 1
    
    def getOrderedClusters(self):
        orderedClusters = []
        cur = self.clusterLLHead
        while cur is not None:
            if len(self.clusters[cur]) > 0:
                orderedClusters.append(list(self.clusters[cur]))
            prv, cur = self.clusterLL[cur]
        #for c in orderedClusters:
        #    print(c)
        return orderedClusters