'''
Created on May 28, 2020

@author: Vlad
'''

import os
import random
import time

from magus_align.decompose import initial_tree, kmh
from magus_helpers import treeutils, sequenceutils
from magus_configuration import Configs

'''
Handles the different ways to decompose the dataset into subsets.
The main way is to estimate a guide tree, then use PASTA's centroid edge decomposition
on the guide tree. Can also decompose randomly (for high speed on huge datasets).
'''

def decomposeSequences(context):
    time1 = time.time()
    
    if len(context.subsetPaths) > 0:
        Configs.log("Subset paths already provided, skipping decomposition..")
    
    elif len(context.subalignmentPaths) > 0:
        context.subsetPaths = context.subalignmentPaths
        Configs.log("Subalignment paths already provided, skipping decomposition..")
    
    else:
        subsetsDir = os.path.join(context.workingDir, "decomposition")
        context.subsetPaths = []
        n = 1
        while True:
            filePath = os.path.join(subsetsDir, "subset_{}.txt".format(n))
            if not os.path.exists(filePath):
                break
            Configs.log("Detected existing subset file {}".format(filePath))
            context.subsetPaths.append(filePath)
            n = n + 1
        
        if len(context.subsetPaths) == 0:
            buildDecomposition(context, subsetsDir)
    
    time2 = time.time()  
    Configs.log("Decomposed {} into {} subsets in {} sec..".format(context.sequencesPath, len(context.subsetPaths), time2-time1))

def buildDecomposition(context, subsetsDir):  
    if not os.path.exists(subsetsDir):
        os.makedirs(subsetsDir)  
    if context.unalignedSequences is None:
        context.unalignedSequences = sequenceutils.readFromFasta(context.sequencesPath, removeDashes=True)    
    
    if (Configs.decompositionStrategy == "random" or context.guideTree == "random") and Configs.outputPath == context.outputFile:
        context.subsetPaths = randomDecomposition(subsetsDir, context.unalignedSequences, Configs.decompositionMaxNumSubsets)
        
    elif Configs.decompositionStrategy == "kmh":
        Configs.log("Decomposing {} with KMH..".format(context.sequencesPath))
        Configs.log("Targetting {} subsets..".format(Configs.decompositionMaxNumSubsets))
        context.subsetPaths = kmh.buildSubsetsKMH(context, subsetsDir)
    
    else:
        guideTreePath  = initial_tree.buildInitialTree(context, subsetsDir, context.guideTree)
        Configs.log("Using target subset size of {}, and maximum number of subsets {}..".format(Configs.decompositionMaxSubsetSize, Configs.decompositionMaxNumSubsets))
        context.subsetPaths = treeutils.decomposeGuideTree(subsetsDir, context.sequencesPath, guideTreePath, 
                                                   Configs.decompositionMaxSubsetSize, Configs.decompositionMaxNumSubsets)        

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
