'''
Created on May 29, 2020

@author: Vlad
'''

import os
import shutil
import time
import random

from magus_align.decompose import decomposer
from magus_helpers import sequenceutils, treeutils, hmmutils
from magus_tasks import task
from magus_tools import external_tools
from magus_configuration import Configs

'''
Not really used for anything, may be removed in the future.
'''

def buildSubsetsKMH(context, subsetsDir):
    tempDir = os.path.join(subsetsDir, "initial_tree")
    
    Configs.log("Building KMH decomposition on {} with skeleton size {}/{}..".format(context.sequencesPath, Configs.decompositionSkeletonSize, 1000))
    time1 = time.time()
    
    initialTreePath, initialAlignPath, unusedTaxa = buildInitialTreeAlign(tempDir, context.sequencesPath) 
    
    if len(unusedTaxa) == 0:
        subsetPaths = treeutils.decomposeGuideTree(tempDir, initialAlignPath, initialTreePath, Configs.decompositionMaxSubsetSize, Configs.decompositionMaxNumSubsets)
    else:
        subsetSeedDir = os.path.join(subsetsDir, "seed_subsets")
        if not os.path.exists(subsetSeedDir):
            os.makedirs(subsetSeedDir)
        subsetSeedPaths = treeutils.decomposeGuideTree(subsetSeedDir, initialAlignPath, initialTreePath, None, Configs.decompositionMaxNumSubsets)
        subsetPaths = reassignTaxons(subsetsDir, subsetSeedPaths, context.unalignedSequences, unusedTaxa)
    
    time2 = time.time()
    Configs.log("Built KMH decomposition on {} in {} sec..".format(context.sequencesPath, time2-time1))

    return subsetPaths

def buildInitialTreeAlign(tempDir, sequencesPath):
    outputTreePath = os.path.join(tempDir, "initial_tree.tre")
    outputAlignPath = os.path.join(tempDir, "initial_align.txt")
    
    if os.path.exists(outputTreePath) and os.path.exists(outputAlignPath):
        return outputTreePath, outputAlignPath
    if os.path.exists(tempDir):
        shutil.rmtree(tempDir)
    os.makedirs(tempDir)
    
    initialAlign, unusedTaxa = decomposer.initial_tree.buildInitialAlignment(sequencesPath, tempDir, Configs.decompositionSkeletonSize, 1000)
    sequenceutils.writeFasta(initialAlign, outputAlignPath)    
    #external_tools.runRaxmlNg(outputAlignPath, tempDir, outputTreePath, 8).run()
    external_tools.runFastTree(outputAlignPath, tempDir, outputTreePath).run()

    return outputTreePath, outputAlignPath, unusedTaxa

def reassignTaxons(subsetsDir, subsetSeedPaths, sequences, unusedTaxa):
    unusedPath = os.path.join(subsetsDir, "unassigned_sequences.txt")
    sequenceutils.writeFasta(sequences, unusedPath, unusedTaxa)
    
    hmmMap = {}
    for subsetPath in subsetSeedPaths:
        hmmDir = os.path.join(os.path.dirname(subsetPath), "hmm_{}".format(os.path.basename(subsetPath)).replace(".", "_"))
        if not os.path.exists(hmmDir):
            os.makedirs(hmmDir)
        hmmMap[subsetPath] = os.path.join(hmmDir, "hmm_model.txt") 
    hmmTasks = hmmutils.buildHmms(hmmMap)
    task.submitTasks(hmmTasks)
    task.awaitTasks(hmmTasks)
    hmmPaths = [t.outputFile for t in hmmTasks]
    
    scoreFileHmmFileMap = {}
    scoreTasks = hmmutils.buildHmmScores(hmmPaths, unusedPath, scoreFileHmmFileMap)
    task.submitTasks(scoreTasks) 
    
    bestScores = {}
    taxonHmmMap = {}
    for scoreTask in task.asCompleted(scoreTasks):
        subsetScores = hmmutils.readSearchFile(scoreTask.outputFile)
        for taxon, scores in subsetScores.items():
            if scores[1] > bestScores.get(taxon, -float("inf")):
                bestScores[taxon] = scores[1]
                taxonHmmMap[taxon] = scoreFileHmmFileMap[scoreTask.outputFile]
    
    subsetTaxons = {file : [] for file in hmmPaths}
    for taxon, hmmPath in taxonHmmMap.items():
        subsetTaxons[hmmPath].append(taxon)
    for subsetPath, hmmPath in hmmMap.items():
        subset = sequenceutils.readFromFasta(subsetPath)
        for taxon in subset:
            subsetTaxons[hmmPath].append(taxon)
    
    subsetPaths = []    
    i = 1
    for hmmPath, subset in subsetTaxons.items():
        subsetPath = os.path.join(subsetsDir, "subset_{}.txt".format(i))
        subsetPaths.append(subsetPath)
        sequenceutils.writeFasta(sequences, subsetPath, subset)
        i = i + 1

    return subsetPaths


    
