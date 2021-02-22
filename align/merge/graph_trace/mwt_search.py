'''
Created on Jun 25, 2020

@author: Vlad
'''

import heapq
from collections import deque 

from configuration import Configs


'''
Resolve clusters into a trace by looking for cycles and removing edges to break the cycles.
We're done when there are no more cycles.
'''

def mwtGreedySearch(graph):
    Configs.log("Finding graph trace with MWT greedy search..")
    
    k = len(graph.context.subalignments)
    lowerBound = [graph.subsetMatrixIdx[i] for i in range(k)]
    upperBound = [graph.subsetMatrixIdx[i] + graph.subalignmentLengths[i] for i in range(k)] 
    
    if graph.clusters is None or len(graph.clusters) == 0:
        graph.buildNodeEdgeDataStructure()
    else:
        graph.buildNodeEdgeDataStructureFromClusters()
    
    context = MwtSearchContext(lowerBound, upperBound)
    state = MwtSearchState()
    state.frontier = list(lowerBound)    
    clusters, totalCost, cycles = greedySearch(graph, state, context)
    graph.clusters = clusters

def mwtSearch(graph):
    Configs.log("Finding graph trace with MWT heuristic search..")
    
    k = len(graph.context.subalignments)
    lowerBound = [graph.subsetMatrixIdx[i] for i in range(k)]
    upperBound = [graph.subsetMatrixIdx[i] + graph.subalignmentLengths[i] for i in range(k)] 
    
    if graph.clusters is None or len(graph.clusters) == 0:
        graph.buildNodeEdgeDataStructure()
    else:
        graph.buildNodeEdgeDataStructureFromClusters()
    
    clusters, totalCost = mwtHeuristicSearch(graph, lowerBound, upperBound)
    graph.clusters = clusters
    
def mwtHeuristicSearch(graph, lowerBound, upperBound):  
    context = MwtSearchContext(lowerBound, upperBound)
    startState = MwtSearchState()
    startState.frontier = list(lowerBound)
    
    heap = []
    visited = set()    
    maxFrontierState = startState
    heapq.heappush(heap, (startState.getHeuristic(), startState))
    
    while len(heap) > 0:
        
        if len(heap) > Configs.searchHeapLimit:
            Configs.log("Heap limit exceeded, clearing heap and moving to max frontier..")
            heap = []
            visited = set()
            heapq.heappush(heap, (maxFrontierState.getHeuristic(), maxFrontierState))
                
        heuristic, state = heapq.heappop(heap)
 
        newMax = False
        newFull = True
        for i in range(len(state.frontier)): 
            if state.frontier[i] > context.maxFrontier[i]:
                newMax = True
            elif state.frontier[i] == context.maxFrontier[i]:
                newFull = False
            elif state.frontier[i] < context.maxFrontier[i]:
                newMax = False
                newFull = False
                break
        if newMax:
            context.maxFrontier = list(state.frontier)
            maxFrontierState = state
            #Configs.log("New frontier {}..".format(graph_partition.cutString(graph, context.maxFrontier)))
            #Configs.log("{} cost, {} removed..".format(state.cost, len(state.removed)))
        if newFull:
            context.fullFrontier = list(state.frontier)
            #Configs.log("New full frontier {}..".format(graph.cutString(context.fullFrontier)))
            heap = []
            visited = set()
         
        percent = getBoundPercent(state.frontier, context)
        if int(percent/10) > int(context.percentDone/10):
            context.percentDone = percent
            Configs.log("{}% done, {} cost, {} edges removed..".format(percent, state.cost, len(state.removed)))
            Configs.log("Max frontier {}..".format(graph.cutString(context.maxFrontier)))
        
        moves = findMoves(graph, state, context, 1)
        if len(moves) == 0:  
            state.frontier = list(lowerBound)          
            isCycle, orderedClusters = findCycleOrCluster(graph, state, context)
            return orderedClusters, state.cost
        
        for nextState in moves:
            stateKey = nextState.getStateKey()
            if stateKey not in visited:
                visited.add(stateKey)
                heapq.heappush(heap, (nextState.getHeuristic(), nextState))        


    Configs.log("Heap empty, resorting to greedy search..")
    context = MwtSearchContext(lowerBound, upperBound)
    state = MwtSearchState()
    state.frontier = list(lowerBound)
    clusters, totalCost, cycles = greedySearch(graph, state, context)
    return clusters, totalCost

def findMoves(graph, state, context, frontierSearchDepth = 3):
    if len(state.frontierRemoved) >= frontierSearchDepth:
        findGreedyProgress(graph, state, context)
    
    oldFrontier = tuple(state.frontier)    
    isCycle, result = findCycleOrCluster(graph, state, context)
                
    if not isCycle:
        return []
    else:
        moves = []
        sameFrontier = oldFrontier == tuple(state.frontier)
        for edge in result:
            newState = MwtSearchState(state)
            newState.removed.add(edge)
            newState.cost = newState.cost + graph.matrix[edge[0]][edge[1]]
            
            if sameFrontier:
                newState.frontierRemoved.add(edge)
            else:
                newState.frontierRemoved = set([edge])
            moves.append(newState)
        return moves

def findGreedyProgress(graph, state, context):
    state.frontierRemoved = set()
    oldFrontier = tuple(state.frontier)
    
    while True:
        isCycle, result = findCycleOrCluster(graph, state, context)

        if oldFrontier != tuple(state.frontier):
            return
            
        if not isCycle:
            return

        edge = min(result, key = lambda x: graph.matrix[x[0]][x[1]])
        state.removed.add(edge)
        state.cost = state.cost + graph.matrix[edge[0]][edge[1]]

def greedySearch(graph, state, context):
    cycles = []
    isCycle, result = findCycleOrCluster(graph, state, context)
    while isCycle:
        cycles.append(result)
        
        percent = getBoundPercent(state.frontier, context)
        if int(percent/10) > int(context.percentDone/10):
            context.maxFrontier = list(state.frontier)
            context.percentDone = percent
            Configs.log("{}% done, {} cost, {} edges removed..".format(percent, state.cost, len(state.removed)))
            Configs.log("Frontier {}..".format(graph.cutString(context.maxFrontier)))
            
        edge = min(result, key = lambda x: graph.matrix[x[0]][x[1]])
        state.removed.add(edge)
        state.cost = state.cost + graph.matrix[edge[0]][edge[1]]
        
        isCycle, result = findCycleOrCluster(graph, state, context)
        
    state.frontier = list(context.lowerBound)
    isCycle, result = findCycleOrCluster(graph, state, context)
    return result, state.cost, cycles
    
def findCycleOrCluster(graph, state, context):    
    k = len(graph.context.subalignments)
    
    curCluster = None
    clusters = []
    orderedClusters = []
    backPointers = {}
    nodeClusters = {}
    
    while True:
        lowerNode, upperNode = None, None
        if curCluster is None:
            for i in range(k):
                if state.frontier[i] < context.upperBound[i]:
                    lowerNode = state.frontier[i]
                    break
            if lowerNode is None:
                return False, orderedClusters
        else:
            for node in clusters[curCluster]:
                asub, apos = graph.matSubPosMap[node]
                if node > state.frontier[asub]:
                    lowerNode = state.frontier[asub]
                    upperNode = node
                    
                    if lowerNode in nodeClusters:
                        #print("Complex cycle..")
                        other = nodeClusters[lowerNode]
                        path = []                        
                        cur = curCluster
                        startNode = upperNode
                        while cur != other:
                            prev, lower, upper = backPointers[cur]
                            segment = findPathBFS(graph, state.frontier, context.upperBound, state.removed, startNode, lower)
                            path.extend(segment)
                            cur = prev
                            startNode = upper
                        segment = findPathBFS(graph, state.frontier, context.upperBound, state.removed, startNode, lowerNode)
                        return True, path
            
            if lowerNode is None:
                for node in clusters[curCluster]:
                    asub, apos = graph.matSubPosMap[node]
                    state.frontier[asub] = state.frontier[asub] + 1
                clusters[curCluster].sort()
                orderedClusters.append(clusters[curCluster])
                prev, lower, upper = backPointers.get(curCluster, (None, 0, 0))
                curCluster = prev
                continue
            
        isCycle, result = findCycleOrClusterFromNode(graph, state, context, lowerNode)
        if isCycle:
            return True, result
        else:
            clusters.append(result)
            idx = len(clusters)-1    
            for node in result:
                nodeClusters[node] = idx
            backPointers[idx] = (curCluster, lowerNode, upperNode) 
            curCluster = idx                                             
    return False, orderedClusters
    
def findCycleOrClusterFromNode(graph, state, context, node):
    k = len(graph.context.subalignments)
    i, pos = graph.matSubPosMap[node]

    stack = [node]
    clusterNodes = set([node])
    seqPositions = {i : node}
    backPointers = {}
    
    while len(stack) > 0:
        curNode = stack.pop()
        for j in range(k):
            sibling = None
            for nbr, value in graph.nodeEdges[curNode][j]:
                if nbr in clusterNodes or edge(curNode, nbr) in state.removed or nbr < state.frontier[j] or nbr >= context.upperBound[j]:
                    continue
                if sibling is not None:
                    return True, [edge(sibling, curNode), edge(curNode, nbr)]
                sibling = nbr
                
                clusterNodes.add(nbr)
                backPointers[nbr] = curNode
                if j in seqPositions:
                    path = set()
                    cur = nbr
                    while cur in backPointers:
                        prev = backPointers[cur]
                        path.add(edge(cur, prev))
                        cur = prev
                    cur = seqPositions[j]
                    while cur in backPointers:
                        prev = backPointers[cur]
                        e = edge(cur, prev)
                        if e in path:
                            path.remove(e)
                        else:
                            path.add(e)
                        cur = prev
        
                    return True, list(path)           
                
                seqPositions[j] = nbr
                stack.append(nbr)  
    
    return False, list(clusterNodes)            

def findPathBFS(graph, lowerBound, upperBound, removed, nodeA, nodeB):
    k = len(graph.context.subalignments)
    asub, apos = graph.matSubPosMap[nodeA]
    bsub, bpos = graph.matSubPosMap[nodeB]
    #queue = deque([nodeA])
    queue = deque([(nodeA, set([asub]))])
    visited = set([nodeA])
    backPointers = {}
    
    while len(queue) > 0:
        curNode, levels = queue.popleft()
        for j in range(k):
            if j in levels and j != bsub:
                continue
            
            for nbr, value in graph.nodeEdges[curNode][j]:
                if nbr in visited or edge(curNode, nbr) in removed or nbr < lowerBound[j] or nbr >= upperBound[j]:
                    continue
                
                visited.add(nbr)
                backPointers[nbr] = curNode                
                #stack.append(nbr)
                
                newlvls = set(levels)
                newlvls.add(j)
                queue.append((nbr, newlvls))
                
                if nbr == nodeB:
                    path = []
                    cur = nbr
                    while cur != nodeA:
                        prev = backPointers[cur]
                        path.append(edge(cur, prev))
                        cur = prev
        
                    return path        
    return None

def edge(a, b):
    return (min(a, b), max(a, b))

def getBoundPercent(bound, context):
    return int(100 * sum([bound[i] - context.lowerBound[i] for i in range(len(bound))]) / context.numNodes)

class MwtSearchContext:
    
    def __init__(self, lowerBound, upperBound):
        self.numNodes = sum([upperBound[i] - lowerBound[i] for i in range(len(lowerBound))])
        self.percentDone = 0
        self.maxFrontier = list(lowerBound)
        self.fullFrontier = list(lowerBound)
        self.lowerBound = list(lowerBound)
        self.upperBound = list(upperBound)       
        
        
class MwtSearchState:
    
    stateCounter = 0
        
    def __init__(self, otherState = None):
        self.frontier = []
        self.removed = set()
        self.frontierRemoved = set()
        self.cost = 0
        
        
        MwtSearchState.stateCounter = MwtSearchState.stateCounter + 1
        self.count = MwtSearchState.stateCounter
        
        if otherState is not None:
            self.frontier = list(otherState.frontier)
            self.removed = set(otherState.removed)
            self.frontierRemoved = set(otherState.frontierRemoved)
            self.cost = otherState.cost
    
       
    def getHeuristic(self):
        return (self.cost, -self.count)
    
    def getHeuristic2(self, lowerBound, upperBound):
        lsum = sum([self.frontier[i] - lowerBound[i] for i in range(len(self.frontier))])
        rsum = sum([upperBound[i] - self.frontier[i] for i in range(len(self.frontier))])
        #return (self.cost/max(lsum,1), -self.count)
        return (self.cost*(1 + rsum/max(lsum,1)), -self.count)
    
    def getStateKey(self):
        #return tuple(self.frontier + [self.cost, len(self.removed)])
        #return tuple(self.frontier)
        return (tuple(self.frontier), frozenset(self.frontierRemoved))
        #return (tuple(self.frontier), self.cost, len(frozenset(self.frontierRemoved)))

'''
def buildEdgeSets(graph):
    k = len(graph.subalignments)
    edgeSets = []
    for node in graph.nodeEdges:
        for i in range(k):
            if len(graph.nodeEdges[node][i]) > 1:
                edgeSets.append([edge(node, nbr) for nbr, value in graph.nodeEdges[node][i]])
                
    return edgeSets

def ilpStage(graph):
    edgeSets = buildEdgeSets(graph)
    removed, cost = edgeSetsAndCyclesToIlp(graph, edgeSets, [])
    Configs.log("Preprocessed with ILP, {} removed, {} cost..".format(len(removed), cost))
    return removed, cost


def mwtIlp(graph, lowerBound, upperBound):
    context = MwtSearchContext(lowerBound, upperBound)    
    state = MwtSearchState()
    state.frontier = list(lowerBound)       
    bestClusters, bestCost, newCycles = greedySearch(graph, state, context)
    cycles = newCycles 
    Configs.log("Initial greedy search: {} cost, {} new cycles found..".format(bestCost, len(newCycles)))
    
    iteration = 0
    while len(newCycles) > 0 and iteration < 10:
        iteration = iteration + 1 
        removed, cost = edgeSetsAndCyclesToIlp(graph, [], cycles)    
        Configs.log("ILP refinement {}: {} cost over {} cycles..".format(iteration, cost, len(cycles)))
            
        context = MwtSearchContext(lowerBound, upperBound)    
        state = MwtSearchState()
        state.frontier = list(lowerBound)   
        state.removed = removed 
        state.cost = cost
        clusters, newCost, newCycles = greedySearch(graph, state, context)
        Configs.log("Greedy search {}: {} cost, {} new cycles found..".format(iteration, newCost, len(newCycles)))
        cycles.extend(newCycles)
        
        if newCost < bestCost:
            bestCost = newCost
            bestClusters = clusters
        else:
            break
    
    return bestClusters, bestCost

def edgeSetsAndCyclesToIlp(graph, edgeSets, cycles):
    from ortools.linear_solver import pywraplp
    solver = pywraplp.Solver('SolveIntegerProblem', pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)
    #solver.EnableOutput()
    objective = solver.Objective()
    
    edgeVars = {}
    for edgeSet in edgeSets:
        edgeConstraint = solver.Constraint(len(edgeSet)-1, len(edgeSet))
        for edge in edgeSet:
            if edge not in edgeVars:            
                edgeVar = solver.IntVar(0, 1, "Edge_{}".format(edge)) 
                edgeVars[edge] = edgeVar
                objective.SetCoefficient(edgeVar, graph.matrix[edge[0]][edge[1]])
            else:
                edgeVar = edgeVars[edge]
            edgeConstraint.SetCoefficient(edgeVar, 1)
    
    for cycle in cycles:
        cycleConstraint = solver.Constraint(1, len(cycle))
        for edge in cycle:
            if edge not in edgeVars:            
                edgeVar = solver.IntVar(0, 1, "Edge_{}".format(edge)) 
                edgeVars[edge] = edgeVar
                objective.SetCoefficient(edgeVar, graph.matrix[edge[0]][edge[1]])
            else:
                edgeVar = edgeVars[edge]
            cycleConstraint.SetCoefficient(edgeVar, 1)
            
    Configs.log("Number of variables: {} ".format(solver.NumVariables()))
    Configs.log("Number of constraints: {}".format(solver.NumConstraints()))
    objective.SetMinimization() 
    status = solver.Solve()
    cost = solver.Objective().Value()
    Configs.log(status)    
    Configs.log("Result: {}".format(cost))
    
    removed = []
    for edge in edgeVars:
        if edgeVars[edge].solution_value() == 1:
            removed.append(edge)
    return set(removed), cost


def edgeSetsAndCyclesToIlpAlternate(graph, edgeSets, cycles):
    from ortools.sat.python import cp_model
    model = cp_model.CpModel()
    solver = cp_model.CpSolver()
    
    edgeVars = {}
    for edgeSet in edgeSets:
        for edge in edgeSet:
            if edge not in edgeVars:
                edgeVars[edge] = model.NewBoolVar("Edge_{}".format(edge))
        model.Add(sum([edgeVars[edge] for edge in edgeSet]) >= len(edgeSet)-1)
    
    for cycle in cycles:
        for edge in cycle:
            if edge not in edgeVars:        
                edgeVars[edge] = model.NewBoolVar("Edge_{}".format(edge))
        model.Add(sum([edgeVars[edge] for edge in cycle]) >= 1)
    
    model.Minimize(sum([graph.matrix[edge[0]][edge[1]]*edgeVars[edge] for edge in edgeVars])) 
    print(model.ModelStats())
    solver.parameters.log_search_progress = True
    solver.parameters.max_time_in_seconds = 600.0
    status = solver.Solve(model)
    cost = solver.ObjectiveValue()
    print(status)    
    print("Result =", cost)  
    
    removed = []
    for edge in edgeVars:
        if solver.Value(edgeVars[edge]) == 1:
            removed.append(edge)
    return set(removed), cost
'''
    