'''
Created on May 28, 2020

@author: Vlad
'''

import os
import random
import time

from align.decompose import pastastyle, kmh
from helpers import treeutils, sequenceutils
from configuration import Configs


def decomposeSequences(task):
    time1 = time.time()
    
    #baseName = os.path.splitext(os.path.basename(task.sequencesPath))[0]
    #subsetsDir = os.path.join(task.workingDir, "decomposition_{}".format(baseName))
    subsetsDir = os.path.join(task.workingDir, "decomposition")
    if not os.path.exists(subsetsDir):
        os.makedirs(subsetsDir)
    
    sequences = sequenceutils.readFromFasta(task.sequencesPath)    
               
    if task.guideTreePath is not None:
        Configs.log("Decomposing {} with user guide tree {}".format(task.sequencesPath, task.guideTreePath))
        Configs.log("Using target subset size of {}, and maximum number of subsets {}..".format(Configs.decompositionMaxSubsetSize, Configs.decompositionMaxNumSubsets))
        task.subsetPaths = treeutils.decomposeGuideTree(subsetsDir, task.sequencesPath, task.guideTreePath, 
                                                   Configs.decompositionMaxSubsetSize, Configs.decompositionMaxNumSubsets)
    
    elif len(sequences) >= 250000:
        task.subsetPaths = randomDecomposition(subsetsDir, sequences, Configs.decompositionMaxNumSubsets)
    
    elif Configs.decompositionStrategy == "pastastyle":
        Configs.log("Decomposing {} with PASTA-style initial tree..".format(task.sequencesPath))
        Configs.log("Using target subset size of {}, and maximum number of subsets {}..".format(Configs.decompositionMaxSubsetSize, Configs.decompositionMaxNumSubsets))
        guideTreePath, initialAlignPath = pastastyle.buildPastaInitialTree(subsetsDir, task.sequencesPath)
        task.subsetPaths = treeutils.decomposeGuideTree(subsetsDir, task.sequencesPath, guideTreePath,
                                                   Configs.decompositionMaxSubsetSize, Configs.decompositionMaxNumSubsets)
        
    elif Configs.decompositionStrategy == "kmh":
        Configs.log("Decomposing {} with KMH..".format(task.sequencesPath))
        Configs.log("Targetting {} subsets..".format(Configs.decompositionMaxNumSubsets))
        task.subsetPaths = kmh.buildSubsetsKMH(subsetsDir, task.sequencesPath)
    
    time2 = time.time()  
    Configs.log("Decomposed {} into {} subsets in {} sec..".format(task.sequencesPath, len(task.subsetPaths), time2-time1))
    

def chooseSkeletonTaxa(sequences, skeletonSize, mode = "fulllength"):
    allTaxa = list(sequences.keys())
    
    if mode == "fulllength":
        seqLengths = [len(sequences[t].seq) for t in sequences]
        
        #topQuartile = numpy.quantile(seqLengths, 0.75)
        seqLengths.sort()
        topQuartile = seqLengths[int(0.75*(len(seqLengths)-1))]
        
        fullLength = []
        notFullLength = []
        for t in allTaxa:
            if abs(len(sequences[t].seq) - topQuartile) < 0.25 * topQuartile:
                fullLength.append(t)
            else:
                notFullLength.append(t) 
        
        random.shuffle(fullLength)
        random.shuffle(notFullLength)     
        allTaxa = fullLength + notFullLength
    else:
        random.shuffle(allTaxa)
        
    skeletonTaxa = allTaxa[:skeletonSize]
    remainingTaxa = allTaxa[skeletonSize:]
    return skeletonTaxa, remainingTaxa

def randomDecomposition(subsetsDir, sequences, numSubsets):
    allTaxa = list(sequences.keys())
    random.shuffle(allTaxa)
    
    taxonSubsets = [allTaxa[i :: numSubsets] for i in range(numSubsets)]
    subsetPaths = []
    for n, subset in enumerate(taxonSubsets):
        subsetPath = os.path.join(subsetsDir, "subset_{}.txt".format(n+1))
        subsetPaths.append(subsetPath)                    
        sequenceutils.writeFasta(sequences, subsetPath, subset) 
    return subsetPaths
