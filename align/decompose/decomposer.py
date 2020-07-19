'''
Created on May 28, 2020

@author: Vlad
'''

import os
import random

from align.decompose import pastastyle, kmh
from helpers import treeutils
from configuration import Configs


def decomposeSequences(workingDir, sequencesPath):
    baseName = os.path.splitext(os.path.basename(sequencesPath))[0]
    subsetsDir = os.path.join(workingDir, "decomposition_{}".format(baseName))
    if not os.path.exists(subsetsDir):
        os.makedirs(subsetsDir)
               
    if Configs.guideTreePath is not None:
        Configs.log("Decomposing {} with user guide tree {}".format(sequencesPath, Configs.guideTreePath))
        Configs.log("Using target subset size of {}, and maximum number of subsets {}..".format(Configs.decompositionMaxSubsetSize, Configs.decompositionMaxNumSubsets))
        subsetPaths = treeutils.decomposeGuideTree(subsetsDir, sequencesPath, Configs.guideTreePath, 
                                                   Configs.decompositionMaxSubsetSize, Configs.decompositionMaxNumSubsets)
    
    elif Configs.decompositionStrategy == "pastastyle":
        Configs.log("Decomposing {} with PASTA-style initial tree..".format(sequencesPath))
        Configs.log("Using target subset size of {}, and maximum number of subsets {}..".format(Configs.decompositionMaxSubsetSize, Configs.decompositionMaxNumSubsets))
        guideTreePath, initialAlignPath = pastastyle.buildPastaInitialTree(subsetsDir, sequencesPath)
        subsetPaths = treeutils.decomposeGuideTree(subsetsDir, sequencesPath, guideTreePath,
                                                   Configs.decompositionMaxSubsetSize, Configs.decompositionMaxNumSubsets)
        
    elif Configs.decompositionStrategy == "kmh":
        Configs.log("Decomposing {} with KMH..".format(sequencesPath))
        Configs.log("Targetting {} subsets..".format(Configs.decompositionMaxNumSubsets))
        subsetPaths = kmh.buildSubsetsKMH(subsetsDir, sequencesPath)
    
    return subsetPaths

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


