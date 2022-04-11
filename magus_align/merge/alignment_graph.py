'''
Created on Apr 14, 2020

@author: Vlad
'''

import os

from magus_helpers import sequenceutils
from magus_configuration import Configs
import threading


'''
Data structure for dealing with alignment graphs.
Subalignment columns are mapped to graph nodes, represented by integers.
Integer nodes can be converted back to corresponding subalignment columns.
Reads/writes graph and cluster files.
'''

class AlignmentGraph:
    
    def __init__(self, context):
        self.context = context
        self.workingDir = os.path.join(self.context.workingDir, "graph")
        self.graphPath = os.path.join(self.workingDir, "graph.txt")
        self.clusterPath = os.path.join(self.workingDir, "clusters.txt")
        self.tracePath = os.path.join(self.workingDir, "trace.txt")
        if not os.path.exists(self.workingDir):
            os.makedirs(self.workingDir)
        
        self.subalignmentLengths = []
        self.subsetMatrixIdx = []
        self.matSubPosMap = []
        
        self.matrixSize = 0
        self.matrix = None
        self.matrixLock = threading.Lock()
        self.nodeEdges = None
        
        self.clusters = []
        self.insertions = set()
        
    def initializeMatrix(self):
        if Configs.constrain:
            self.subalignmentLengths = [sequenceutils.readSequenceLengthFromFasta(file) for file in self.context.subalignmentPaths]
        else:
            self.subalignmentLengths = [len(self.context.unalignedSequences[s[0]].seq) for s in self.context.subalignments]
        
        self.matrixSize = sum(self.subalignmentLengths)    
        self.subsetMatrixIdx = [0] * len(self.subalignmentLengths)
        for k in range(1, len(self.subalignmentLengths)):        
            self.subsetMatrixIdx[k] = self.subsetMatrixIdx[k-1] + self.subalignmentLengths[k-1]
        
        self.matSubPosMap = [0] * self.matrixSize
        i = 0
        for k in range(len(self.subalignmentLengths)):
            for j in range(self.subalignmentLengths[k]):
                self.matSubPosMap[i] = (k, j)
                i = i + 1
        
        self.matrix = [{} for i in range(self.matrixSize)]
    
    def writeGraphToFile(self, filePath):
        with open(filePath, 'w') as textFile:
            for i in range(len(self.matrix)):
                for k in self.matrix[i]:
                    textFile.write("{} {} {}\n".format(i, k, self.matrix[i][k]))
        Configs.log("Wrote matrix to {}".format(filePath))

    def readGraphFromFile(self, filePath):
        self.matrix = [{} for i in range(self.matrixSize)]
        with open(filePath) as f:
            for line in f:
                tokens = [int(token) for token in line.strip().split()]
                self.matrix[tokens[0]][tokens[1]] = tokens[2]
        Configs.log("Read matrix from {}".format(filePath))
    
    def writeClustersToFile(self, filePath):
        with open(filePath, 'w') as textFile:
            for cluster in self.clusters:
                textFile.write("{}\n".format(" ".join([str(c) for c in cluster])))
        
    def readClustersFromFile(self, filePath):
        self.clusters = []
        with open(filePath) as f:
            for line in f:
                tokens = [int(token) for token in line.strip().split()]
                if len(tokens) > 1:
                    self.clusters.append(tokens) 
        print("Found {} clusters..".format(len(self.clusters)))
    
    def buildNodeEdgeDataStructure(self):
        Configs.log("Preparing node edge data structure..")
        k = len(self.subalignmentLengths)
        self.nodeEdges = {}
        
        for a in range(self.matrixSize):
            asub, apos = self.matSubPosMap[a]
            self.nodeEdges[a] = [[] for i in range(k)]
            for b, value in self.matrix[a].items():
                bsub, bpos = self.matSubPosMap[b] 
                if asub == bsub:
                    continue
                self.nodeEdges[a][bsub].append((b, value))
            for i in range(k):
                self.nodeEdges[a][i].sort(key = lambda pair: pair[0])
        Configs.log("Prepared node edge data structure..")
    
    def buildNodeEdgeDataStructureFromClusters(self):
        Configs.log("Preparing node edge data structure..")
        k = len(self.subalignmentLengths)
        self.nodeEdges = {}
        
        Configs.log("Using {} pre-existing clusters to simplify alignment graph..".format(len(self.clusters)))
        for a in range(self.matrixSize):
            self.nodeEdges[a] = [[] for i in range(k)]
        
        for cluster in self.clusters:
            for a in cluster:
                asub, apos = self.matSubPosMap[a]
                for b in cluster:
                    bsub, bpos = self.matSubPosMap[b]
                    if asub == bsub or b not in self.matrix[a]:
                        continue
                    value = self.matrix[a][b]
                    self.nodeEdges[a][bsub].append((b, value))
                for i in range(k):
                    self.nodeEdges[a][i].sort(key = lambda pair: pair[0])
        Configs.log("Prepared node edge data structure..")

    def cutString(self, cut):
        stringCut = list(cut)
        for i, value in enumerate(stringCut):
            stringCut[i] = value - self.subsetMatrixIdx[i]
        return stringCut

    def computeClusteringCost(self, clusters):
        cutCost = 0
        nodeClusters = {}
        
        for n, cluster in enumerate(clusters):
            for a in cluster:
                nodeClusters[a] = n
                
        clusterCounter = len(clusters)
        for a in range(self.matrixSize):
            if a not in nodeClusters:
                nodeClusters[a] = clusterCounter
                clusterCounter = clusterCounter + 1
                 
        for a in range(self.matrixSize):
            asub, apos = self.matSubPosMap[a]
            for b, value in self.matrix[a].items():
                bsub, bpos = self.matSubPosMap[b] 
                if asub != bsub and nodeClusters[a] != nodeClusters[b]:
                    cutCost = cutCost + value
    
        return int(cutCost/2) 
    
    def addSingletonClusters(self):
        newClusters = []
        
        lastIdx = list(self.subsetMatrixIdx) 
        for cluster in self.clusters:
            for a in cluster:
                #print(a)
                asub, apos = self.matSubPosMap[a]                
                for node in range(lastIdx[asub], a):
                    newClusters.append([node])
                lastIdx[asub] = a+1
            newClusters.append(cluster)
        for i in range(len(lastIdx)):
            for node in range(lastIdx[i], self.subsetMatrixIdx[i] + self.subalignmentLengths[i]):
                newClusters.append([node])
        self.clusters = newClusters
        return newClusters
    
    
    
    