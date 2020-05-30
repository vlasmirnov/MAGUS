'''
Created on May 29, 2020

@author: Vlad
'''

import os
import shutil
import time

from align.decompose import decomposer
from helpers import sequenceutils, treeutils, hmmutils
from tools import external_tools
from configuration import Configs


def buildSubsetsKMH(subsetsDir, sequencesPath):
    tempDir = os.path.join(subsetsDir, "initial_tree")
    
    Configs.log("Building KMH decomposition on {} with skeleton size {}..".format(sequencesPath, Configs.decompositionSkeletonSize))
    time1 = time.time()
    
    initialTreePath, initialAlignPath = buildInitialTreeAlign(tempDir, sequencesPath) 
    subsetSeedPaths = treeutils.decomposeGuideTree(tempDir, initialAlignPath, initialTreePath, None, Configs.decompositionMaxNumSubsets)
    subsetPaths = reassignTaxons(subsetsDir, subsetSeedPaths, sequencesPath)
    
    time2 = time.time()
    Configs.log("Built KMH decomposition on {} in {} sec..".format(sequencesPath, time2-time1))

    return subsetPaths

def buildInitialTreeAlign(tempDir, sequencesPath):
    outputTreePath = os.path.join(tempDir, "initial_tree.tre")
    outputAlignPath = os.path.join(tempDir, "initial_align.txt")
    
    if os.path.exists(outputTreePath) and os.path.exists(outputAlignPath):
        return outputTreePath, outputAlignPath
    if os.path.exists(tempDir):
        shutil.rmtree(tempDir)
    os.makedirs(tempDir)
    
    skeletonPath = os.path.join(tempDir, "skeleton_sequences.txt")
    
    sequences = sequenceutils.readFromFasta(sequencesPath, True)
    skeletonTaxa, remainingTaxa = decomposer.chooseSkeletonTaxa(sequences, Configs.decompositionSkeletonSize)
    sequenceutils.writeFasta(sequences, skeletonPath, skeletonTaxa)
    external_tools.runMafft(skeletonPath, None, tempDir, outputAlignPath, Configs.numCores)
    
    external_tools.runRaxmlNg(outputAlignPath, tempDir, outputTreePath, 8)
    #external_tools.runFastTree(outputAlignPath, tempDir, outputTreePath)

    return outputTreePath, outputAlignPath

def reassignTaxons(subsetsDir, subsetSeedPaths, sequencesPath):
    sequences = sequenceutils.readFromFasta(sequencesPath, True)
    hmmPaths = {}
    for subsetPath in subsetSeedPaths:
        hmmPaths[subsetPath] = os.path.join(os.path.dirname(subsetPath), "hmm_{}".format(os.path.basename(subsetPath)).replace(".", "_"))
    hmmutils.buildHmms(hmmPaths)
    subsetScores = hmmutils.buildHmmScores(hmmPaths, sequencesPath)
    
    #for subsetPath in hmmPaths:
    #    shutil.rmtree(hmmPaths[subsetPath])
    
    bestScores = {}
    taxonHmmMap = {}
    for i, subsetPath in enumerate(subsetSeedPaths):
        scores = subsetScores[subsetPath]
        for taxon in scores:
            if scores[taxon] > bestScores.get(taxon, -float("inf")):
                bestScores[taxon] = scores[taxon]
                taxonHmmMap[taxon] = i
    
    subsetTaxons = [[] for subsetPath in subsetSeedPaths]
    for taxon in taxonHmmMap:
        subsetTaxons[taxonHmmMap[taxon]].append(taxon)
    
    subsetPaths = []    
    for i, subset in enumerate(subsetTaxons):
        subsetPath = os.path.join(subsetsDir, "subset_{}.txt".format(i+1))
        subsetPaths.append(subsetPath)
        sequenceutils.writeFasta(sequences, subsetPath, subset)

    return subsetPaths


    
