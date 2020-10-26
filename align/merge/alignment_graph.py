'''
Created on Apr 14, 2020

@author: Vlad
'''

import os

from helpers import sequenceutils, tasks
from configuration import Configs

class AlignmentGraph:
    
    
    def __init__(self, task):
        self.workingDir = os.path.join(task.workingDir, "graph")
        self.graphPath = os.path.join(self.workingDir, "graph.txt")
        if not os.path.exists(self.workingDir):
            os.makedirs(self.workingDir)
        
        self.subalignmentPaths = None
        self.clusterPath = None
        
        self.unalignment = {}
        self.subalignments = []
        self.taxonSubalignmentMap = {}
        self.subalignmentLengths = []
        self.subsetMatrixIdx = []
        self.matSubPosMap = []
        
        self.matrixSize = 0
        self.matrix = None
        self.matrixLock = None
        self.nodeEdges = None
        
        self.clusters = []
        
        
    def loadSubalignments(self, subalignmentPaths):
        self.subalignmentPaths = subalignmentPaths
        self.subalignments = [sequenceutils.readFromFasta(filePath) for filePath in subalignmentPaths]
        for i, subalignment in enumerate(self.subalignments):              
            for taxon in subalignment:
                self.taxonSubalignmentMap[taxon] = i
                self.unalignment[taxon] = sequenceutils.Sequence(taxon, subalignment[taxon].seq.replace('-',''))
        
        self.subalignmentLengths = [len(next(iter(subset.values())).seq) for subset in self.subalignments]
        self.subsetMatrixIdx = [0] * len(self.subalignments)
        for k in range(1, len(self.subalignments)):        
            self.subsetMatrixIdx[k] = self.subsetMatrixIdx[k-1] + self.subalignmentLengths[k-1]
        
        self.matrixSize = sum(self.subalignmentLengths)
        
        self.matSubPosMap = [0] * self.matrixSize
        i = 0
        for k in range(len(self.subalignments)):
            for j in range(self.subalignmentLengths[k]):
                self.matSubPosMap[i] = (k, j)
                i = i + 1
    
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
        k = len(self.subalignments)
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
        k = len(self.subalignments)
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
        
    def clustersToAlignment(self, filePath):
        finalAlignment = {}
        for subset in self.subalignments:
            for taxon in subset:
                finalAlignment[taxon] = sequenceutils.Sequence(taxon, ['-'] * len(self.clusters))
        
        for idx, cluster in enumerate(self.clusters):
            for b in cluster:
                bsub, bpos = self.matSubPosMap[b]
                for taxon in self.subalignments[bsub]:
                    finalAlignment[taxon].seq[idx] = self.subalignments[bsub][taxon].seq[bpos]

        for taxon in finalAlignment:
            finalAlignment[taxon].seq = "".join(finalAlignment[taxon].seq)    
                                
        sequenceutils.writeFasta(finalAlignment, filePath)  
        Configs.log("Wrote final output alignment to {}".format(filePath))

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
    
    def cheaterAlignment(self, filePath):    
        for s, subset in enumerate(self.subalignments):
            subsetAlignPath = os.path.join(self.workingDir, os.path.basename(self.subalignmentPaths[s]))
            locker = os.path.join(self.workingDir, "lock_{}".format(os.path.basename(subsetAlignPath)))
            try:
                with open(locker, 'x') as subfile:
                    pass
            except:
                continue
            
            if os.path.exists(subsetAlignPath):
                finalAlignment = sequenceutils.readFromFasta(subsetAlignPath, removeDashes=False)
            else:
                Configs.log("Parsing subset {}..".format(os.path.basename(self.subalignmentPaths[s])))
                finalAlignment = {}
                for taxon in subset:
                    finalAlignment[taxon] = sequenceutils.Sequence(taxon, [])
                
                lastIdx = [0 for i in self.subsetMatrixIdx]
                curLen = [0 for i in self.subsetMatrixIdx]
                for idx, cluster in enumerate(self.clusters):
                    
                    newLen = 0
                    for b in cluster:
                        bsub, bpos = self.matSubPosMap[b]
                        if curLen[bsub] + bpos - lastIdx[bsub] > newLen:
                            newLen = curLen[bsub] + bpos - lastIdx[bsub]
                            
                    for b in cluster:   
                        bsub, bpos = self.matSubPosMap[b]     
                        if bsub == s:
                            for taxon in self.subalignments[bsub]:
                                for n in range(lastIdx[bsub], bpos):
                                    finalAlignment[taxon].seq.append(self.subalignments[bsub][taxon].seq[n])
                                for n in range(newLen - (curLen[bsub] + bpos - lastIdx[bsub])):
                                    finalAlignment[taxon].seq.append('-')
                                finalAlignment[taxon].seq.append(self.subalignments[bsub][taxon].seq[bpos])
                        lastIdx[bsub] = bpos+1
                        curLen[bsub] = newLen+1
                
                newLen = 0
                for bsub, bpos in enumerate(self.subalignmentLengths): 
                    if curLen[bsub] + bpos - lastIdx[bsub]  > newLen:
                        newLen = curLen[bsub] + bpos - lastIdx[bsub]
                
                for taxon in subset:
                    for n in range(lastIdx[s], self.subalignmentLengths[s]):
                        finalAlignment[taxon].seq.append(self.subalignments[s][taxon].seq[n])
                    for n in range(newLen - (curLen[s] + self.subalignmentLengths[s] - lastIdx[s])):
                        finalAlignment[taxon].seq.append('-')
                        
                for taxon in finalAlignment:
                    finalAlignment[taxon].seq = "".join(finalAlignment[taxon].seq) 
                sequenceutils.writeFasta(finalAlignment, subsetAlignPath)
            
                  
            tasks.lockFiles(os.path.join(self.workingDir, "output.lock")) 
            Configs.log("Writing in alignment of length {}".format(len(next(iter(finalAlignment.values())).seq)))
            sequenceutils.writeFasta(finalAlignment, filePath, append = True) 
            tasks.unlockFiles(os.path.join(self.workingDir, "output.lock"))
                
             
        Configs.log("Wrote final output alignment to {}".format(filePath))
    
    '''
    def cheaterAlignment(self, filePath):    
        finalAlignment = {}
        for subset in self.subalignments:
            for taxon in subset:
                finalAlignment[taxon] = sequenceutils.Sequence(taxon, [])
        
        lastIdx = [0 for i in self.subsetMatrixIdx]
        for idx, cluster in enumerate(self.clusters):
            Configs.log("Parsing cluster {}..".format(idx))
            maxDist = 0
            for b in cluster:
                bsub, bpos = self.matSubPosMap[b]
                if bpos - lastIdx[bsub] > maxDist:
                    maxDist = bpos - lastIdx[bsub]
                    
            for b in cluster:   
                bsub, bpos = self.matSubPosMap[b]     
                for taxon in self.subalignments[bsub]:
                    for n in range(lastIdx[bsub], bpos):
                        finalAlignment[taxon].seq.append(self.subalignments[bsub][taxon].seq[n])
                    for n in range(maxDist - bpos + lastIdx[bsub]):
                        finalAlignment[taxon].seq.append('-')
                    finalAlignment[taxon].seq.append(self.subalignments[bsub][taxon].seq[bpos])
                lastIdx[bsub] = bpos+1
        
        maxDist = 0
        for bsub, bpos in enumerate(self.subalignmentLengths): 
            if bpos - lastIdx[bsub] > maxDist:
                maxDist = bpos - lastIdx[bsub]
        
        for bsub, bpos in enumerate(self.subalignmentLengths):   
            for taxon in self.subalignments[bsub]:
                for n in range(lastIdx[bsub], bpos):
                    finalAlignment[taxon].seq.append(self.subalignments[bsub][taxon].seq[n])
                for n in range(maxDist - bpos + lastIdx[bsub]):
                    finalAlignment[taxon].seq.append('-')


        for taxon in finalAlignment:
            finalAlignment[taxon].seq = "".join(finalAlignment[taxon].seq)    
                                
        sequenceutils.writeFasta(finalAlignment, filePath)  
        Configs.log("Wrote final output alignment to {}".format(filePath))
    '''    
    '''
    def plotClusterHistogram(self, clusters, numBins, outputFile):
        import matplotlib.pyplot as plt
        from matplotlib import colors
        fig, ax = plt.subplots(1, 1, tight_layout=True)
        plt.rcParams.update({'font.size': 15})
        
        totalClusters = clusters[:]
        nodeClusters = {}
        for n, cluster in enumerate(clusters):
            for a in cluster:
                nodeClusters[a] = n
        for a in range(self.matrixSize):
            if a not in nodeClusters:
                totalClusters.append([a])
        print("Plotting out {} clusters..".format(len(totalClusters)))
        print("Recomputed clustering cost {}".format(self.computeClusteringCost(totalClusters)))
        x = [len(c) for c in totalClusters]
        
        ax.set_facecolor((0.9,0.9,0.9))        
        N, bins, patches = plt.hist(x, bins=numBins)

        fracs = N / N.max()        
        norm = colors.Normalize(fracs.min(), fracs.max())
        
        for thisfrac, thispatch in zip(fracs, patches):
            color = plt.cm.magma(norm(thisfrac))
            #thispatch.set_facecolor(color)
        
        fig.canvas.draw()
        #fig.savefig(outputFile, bbox_inches = 'tight', pad_inches = 0, dpi=300)
        fig.savefig(outputFile, bbox_inches = 'tight', pad_inches = 0)
    '''