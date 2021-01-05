'''
Created on Apr 14, 2020

@author: Vlad
'''

import os

from helpers import sequenceutils
from configuration import Configs
from tasks import task
import threading
import shutil

class AlignmentGraph:
    
    
    def __init__(self, context):
        self.context = context
        self.workingDir = os.path.join(self.context.workingDir, "graph")
        self.graphPath = os.path.join(self.workingDir, "graph.txt")
        self.clusterPath = os.path.join(self.workingDir, "clusters.txt")
        self.tracePath = os.path.join(self.workingDir, "trace.txt")
        self.packedAlignmentPath = os.path.join(self.workingDir, "packed_alignment.txt")
        self.alignmentMaskPath = os.path.join(self.workingDir, "alignment_mask.txt")
        self.maskedSequencesPath = os.path.join(self.workingDir, "masked_sequences.txt")
        if not os.path.exists(self.workingDir):
            os.makedirs(self.workingDir)
        
        self.subalignmentLengths = []
        self.subsetMatrixIdx = []
        self.matSubPosMap = []
        
        self.backbonePaths = []
        self.backboneTasks = []
        self.backboneExtend = set()
        self.backboneTaxa = {}
        self.backboneSubsetAlignment = {}
        
        self.matrixSize = 0
        self.matrix = None
        self.matrixLock = threading.Lock()
        self.nodeEdges = None
        
        self.clusters = []
        
    def initializeMatrix(self):
        self.subalignmentLengths = [sequenceutils.readSequenceLengthFromFasta(file) for file in self.context.subalignmentPaths]
        self.subsetMatrixIdx = [0] * len(self.subalignmentLengths)
        for k in range(1, len(self.subalignmentLengths)):        
            self.subsetMatrixIdx[k] = self.subsetMatrixIdx[k-1] + self.subalignmentLengths[k-1]
        
        self.matrixSize = sum(self.subalignmentLengths)
        
        self.matSubPosMap = [0] * self.matrixSize
        i = 0
        for k in range(len(self.subalignmentLengths)):
            for j in range(self.subalignmentLengths[k]):
                self.matSubPosMap[i] = (k, j)
                i = i + 1
        
        self.matrix = [{} for i in range(self.matrixSize)]
    
    def initializeBackboneSequenceMapping(self):
        if len(self.backboneTaxa) == 0:
            backboneSubsetTaxonMap = self.context.subsetTaxonMap
        else:
            backboneSubsetTaxonMap = {}
            for taxon in self.backboneTaxa:
                i = self.context.taxonSubsetMap[taxon]
                backboneSubsetTaxonMap[i] = backboneSubsetTaxonMap.get(i, [])
                backboneSubsetTaxonMap[i].append(taxon) 
        
        for i, subalignPath in enumerate(self.context.subalignmentPaths):
            subalignment = sequenceutils.readFromFasta(subalignPath, removeDashes=False)
            for taxon in backboneSubsetTaxonMap.get(i, []):
                self.backboneSubsetAlignment[taxon] = subalignment[taxon]
    
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
        k = len(self.context.subalignmentPaths)
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
        k = len(self.context.subalignmentPaths)
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
    
    def clustersToPackedAlignment(self):
        self.addSingletonClusters()
        packedAlignment = {}
        for path in self.context.subalignmentPaths:            
            packedAlignment[path] = sequenceutils.Sequence(path, ['-'] * len(self.clusters))

        for idx, cluster in enumerate(self.clusters):
            for b in cluster:
                bsub, bpos = self.matSubPosMap[b]
                packedAlignment[self.context.subalignmentPaths[bsub]].seq[idx] = 'X'

        for s in packedAlignment:
            packedAlignment[s].seq = "".join(packedAlignment[s].seq)   
                                
        sequenceutils.writeFasta(packedAlignment, self.packedAlignmentPath)  
        Configs.log("Wrote uncompressed packed subset alignment to {}".format(self.packedAlignmentPath))
        if not os.path.exists(self.alignmentMaskPath):  
            self.buildAlignmentMask()
    
    def buildAlignmentMask(self):
        numCells = len(self.context.unalignedSequences) * len(self.clusters)
        Configs.log("Final alignment will have {} cells..".format(numCells))
        Configs.log("Building alignment mask..")
        if numCells > 5e10:
            Configs.log("Removing columns with more than 99% gaps..")
            numSeq = len(self.context.unalignedSequences)
            portion = 0.99
            mask = [0] * len(self.clusters) 
            gapCounts = self.buildGapCounts()
            
            for n, cluster in enumerate(self.clusters):
                nonGaps = 0
                for b in cluster:
                    bsub, bpos = self.matSubPosMap[b]
                    nonGaps = nonGaps + len(self.context.subsetTaxonMap[bsub]) - gapCounts[self.context.subalignmentPaths[bsub]][bpos]
                    if nonGaps >= (1-portion) * numSeq:
                        mask[n] = 1
                        break       
        else:
            mask = [1] * len(self.clusters)
            
        maskString = "".join([str(c) for c in mask])
        with open(self.alignmentMaskPath, 'w') as u:
            u.write("{}\n".format(maskString))
        
    def buildGapCounts(self):
        countingTasks = []
        for subalignPath in self.context.subalignmentPaths:
            countsPath = os.path.join(self.workingDir, "gap_counts_{}".format(os.path.basename(subalignPath)))
            args = {"alignFile" : subalignPath, "outputFile" : countsPath}
            countingTask = task.Task(taskType = "recordGapCounts", outputFile = args["outputFile"], taskArgs = args)
            countingTasks.append(countingTask)
        
        if len(self.context.unalignedSequences) >= 10000:
            task.submitTasks(countingTasks)
        else:
            for countingTask in countingTasks:
                countingTask.run()
                countingTask.isFinished = True
        
        gapCounts = {}
        for countingTask in task.asCompleted(countingTasks):
            with open(countingTask.outputFile) as f:
                countLine = [int(c) for c in f.readline().strip().split()]
                gapCounts[countingTask.taskArgs["alignFile"]] = countLine
            os.remove(countingTask.outputFile)    
        return gapCounts
    
    def writeUnpackedAlignment(self, filePath):
        tempFile = os.path.join(os.path.dirname(filePath), "temp_{}".format(os.path.basename(filePath)))
        if os.path.exists(tempFile):
            os.remove(tempFile)
            
        tempMaskFile = os.path.join(self.workingDir, "temp_{}".format(os.path.basename(self.maskedSequencesPath)))
        if os.path.exists(tempMaskFile):
            os.remove(tempMaskFile)
            
        Configs.log("Assembling final alignment in {}".format(filePath))
        inducedSubalignTasks = []
        for subalignPath in self.context.subalignmentPaths:
            inducedAlignPath = os.path.join(self.workingDir, "induced_{}".format(os.path.basename(subalignPath)))
            maskedSeqPath = os.path.join(self.workingDir, "masked_{}".format(os.path.basename(subalignPath)))
            args = {"packedFilePath" : self.packedAlignmentPath, "subalignmentPath" : subalignPath, "outputFile" : inducedAlignPath, 
                    "alignmentMaskPath" : self.alignmentMaskPath, "maskedSequencesPath" : maskedSeqPath}
            inducedTask = task.Task(taskType = "buildInducedSubalignment", outputFile = args["outputFile"], taskArgs = args)
            inducedSubalignTasks.append(inducedTask)
        
        if len(self.context.unalignedSequences) >= 10000:
            task.submitTasks(inducedSubalignTasks)
        else:
            for inducedTask in inducedSubalignTasks:
                inducedTask.run()
                inducedTask.isFinished = True
                
        for inducedTask in task.asCompleted(inducedSubalignTasks):
            inducedAlign = sequenceutils.readFromFasta(inducedTask.outputFile, removeDashes=False)
            Configs.log("Appending induced alignment, {} sequences of length {}..".format(len(inducedAlign), len(next(iter(inducedAlign.values())).seq)))
            sequenceutils.writeFasta(inducedAlign, tempFile, append = True)   
            
            maskedSeq = sequenceutils.readFromFasta(inducedTask.taskArgs["maskedSequencesPath"])
            sequenceutils.writeFasta(maskedSeq, tempMaskFile, append = True)
              
            os.remove(inducedTask.outputFile)
            os.remove(inducedTask.taskArgs["maskedSequencesPath"])
        shutil.move(tempMaskFile, self.maskedSequencesPath)
        shutil.move(tempFile, filePath)
        Configs.log("Wrote final alignment to {}".format(filePath))    
    
def recordGapCounts(**kwargs):
    counts = sequenceutils.countGaps(kwargs["alignFile"])
    countString = " ".join([str(c) for c in counts])
    with open(kwargs["outputFile"], 'w') as u:
        u.write("{}\n".format(countString))
    
def buildInducedSubalignment(**kwargs):
    packedFilePath = kwargs["packedFilePath"]
    subalignmentPath = kwargs["subalignmentPath"]
    alignmentMaskPath = kwargs["alignmentMaskPath"]
    maskedSequencesPath = kwargs["maskedSequencesPath"]
    inducedAlignPath = kwargs["outputFile"]
    
    tempMaskedSequencesPath = os.path.join(os.path.dirname(maskedSequencesPath), "temp_{}".format(os.path.basename(maskedSequencesPath)))
    tempInducedAlignPath = os.path.join(os.path.dirname(inducedAlignPath), "temp_{}".format(os.path.basename(inducedAlignPath)))
    
    packedAlign = sequenceutils.readFromFasta(packedFilePath, removeDashes=False)
    subsetAlign = sequenceutils.readFromFasta(subalignmentPath, removeDashes=False)
    mask = packedAlign[subalignmentPath].seq
    inducedAlign = {taxon : [] for taxon in subsetAlign}
    maskedAlignment = {taxon : [] for taxon in subsetAlign}
    with open(alignmentMaskPath) as f:
        inclusionMask = [int(c) for c in f.readline().strip()]
        
    j = 0
    for i in range(len(mask)):
        if mask[i] == '-' and inclusionMask[i] == 1:
            for taxon in subsetAlign:
                inducedAlign[taxon].append('-')
        elif mask[i] != '-':
            for taxon in subsetAlign:
                letter = subsetAlign[taxon].seq[j]
                if inclusionMask[i] == 1:
                    inducedAlign[taxon].append(letter)
                if letter != '-':
                    maskedAlignment[taxon].append(letter.upper() if inclusionMask[i] == 1 else letter.lower())
            j = j + 1
    
    for s in maskedAlignment:
        maskedAlignment[s] = sequenceutils.Sequence(s, "".join(maskedAlignment[s]))   
    sequenceutils.writeFasta(maskedAlignment, tempMaskedSequencesPath)
    shutil.move(tempMaskedSequencesPath, maskedSequencesPath)
            
    for s in inducedAlign:
        inducedAlign[s] = sequenceutils.Sequence(s, "".join(inducedAlign[s]))
    sequenceutils.writeFasta(inducedAlign, tempInducedAlignPath)
    shutil.move(tempInducedAlignPath, inducedAlignPath)
    
    