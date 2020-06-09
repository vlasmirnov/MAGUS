'''
Created on Apr 14, 2020

@author: Vlad
'''


from helpers import sequenceutils
from configuration import Configs

class LibraryGraph:
    
    
    def __init__(self, workingDir):
        self.workingDir = workingDir
        self.subsetPaths = None
        self.graphPath = None
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
        
        self.clusters = []
        
        
    def loadSubalignments(self, subsetPaths):
        '''
        self.subalignments = []
        for filePath in subsetPaths:
            baseName = os.path.basename(filePath)
            cleanFilePath = os.path.join(Configs.workingDir, "clean_{}".format(baseName))
            sequenceutils.cleanGapColumns(filePath, cleanFilePath)
            self.subalignments.append(sequenceutils.readFromFasta(cleanFilePath))
        '''
        self.subsetPaths = subsetPaths
        self.subalignments = [sequenceutils.readFromFasta(filePath) for filePath in subsetPaths]
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
        
    def clustersToAlignment(self, filePath):
        finalAlignment = {}
        for subset in self.subalignments:
            for taxon in subset:
                finalAlignment[taxon] = sequenceutils.Sequence(taxon, "")
        
        lastIdx = {}
        for cluster in self.clusters:
            for b in cluster:
                bsub, bpos = self.matSubPosMap[b]
                bposlast = lastIdx.get(bsub, -1)
                
                for p in range(bposlast + 1, bpos):
                    for i in range(len(self.subalignments)):
                        for taxon in self.subalignments[i]:
                            if i == bsub:
                                finalAlignment[taxon].seq = finalAlignment[taxon].seq + self.subalignments[bsub][taxon].seq[p]
                            else:
                                finalAlignment[taxon].seq = finalAlignment[taxon].seq + '-'
                
                lastIdx[bsub] = bpos
            
            clusterSets = set()    
            for b in cluster:
                bsub, bpos = self.matSubPosMap[b]
                clusterSets.add(bsub)
                for taxon in self.subalignments[bsub]:
                    finalAlignment[taxon].seq = finalAlignment[taxon].seq + self.subalignments[bsub][taxon].seq[bpos]
                #lastIdx[bsub] = lastIdx[bsub] + 1
            
            for i in range(len(self.subalignments)):
                if i not in clusterSets:
                    for taxon in self.subalignments[i]:
                        finalAlignment[taxon].seq = finalAlignment[taxon].seq + '-'
                
        
        for bsub in range(len(self.subalignments)):
            bposlast = lastIdx.get(bsub, -1)
                        
            for p in range(bposlast + 1, self.subalignmentLengths[bsub]):
                for i in range(len(self.subalignments)):
                    for taxon in self.subalignments[i]:
                        if i == bsub:
                            finalAlignment[taxon].seq = finalAlignment[taxon].seq + self.subalignments[bsub][taxon].seq[p]
                        else:
                            finalAlignment[taxon].seq = finalAlignment[taxon].seq + '-'        
                                
        sequenceutils.writeFasta(finalAlignment, filePath)  
        Configs.log("Wrote final output alignment to {}".format(filePath))
