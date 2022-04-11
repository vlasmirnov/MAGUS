'''
Created on May 29, 2020

@author: Vlad
'''

import os
import shutil
import time
import random

from magus_align.decompose import decomposer
from magus_helpers import sequenceutils, hmmutils, treeutils
from magus_tasks import task
from magus_tools import external_tools
from magus_configuration import Configs

'''
Different options for estimating a guide tree.
The main ones are FastTree (for accuracy) and Clustal Omega's mbed (for speed).
'''

def buildInitialTree(context, workingDir, treeType):
    if treeType is not None and os.path.exists(treeType):
        Configs.log("Found user guide tree {}".format(treeType))
        return treeType
    
    tempDir = os.path.join(workingDir, "initial_tree")
    outputTreePath = os.path.join(tempDir, "initial_tree.tre")
    if os.path.exists(outputTreePath):
        Configs.log("Found existing initial tree {}".format(outputTreePath))
        return outputTreePath
    if os.path.exists(tempDir):
        shutil.rmtree(tempDir)
    os.makedirs(tempDir)
    
    time1 = time.time() 
    
    if treeType is None or treeType.lower() == "fasttree": 
        Configs.log("Building PASTA-style FastTree initial tree on {} with skeleton size {}..".format(context.sequencesPath, Configs.decompositionSkeletonSize))
        alignPath = os.path.join(tempDir, "initial_align.txt")
        buildInitialAlignment(context.unalignedSequences, tempDir, Configs.decompositionSkeletonSize, None, alignPath)
        external_tools.runFastTree(alignPath, tempDir, outputTreePath, "fast").run()
    elif treeType is None or treeType.lower() == "fasttree-noml": 
        Configs.log("Building PASTA-style FastTree (NO ML) initial tree on {} with skeleton size {}..".format(context.sequencesPath, Configs.decompositionSkeletonSize))
        alignPath = os.path.join(tempDir, "initial_align.txt")
        buildInitialAlignment(context.unalignedSequences, tempDir, Configs.decompositionSkeletonSize, None, alignPath)
        external_tools.runFastTree(alignPath, tempDir, outputTreePath, "noml").run()
    elif treeType.lower() == "parttree":
        Configs.log("Building MAFFT PartTree initial tree on {}..".format(context.sequencesPath))
        taxa = list(context.unalignedSequences.keys())
        external_tools.runMafftGuideTree(context.sequencesPath, tempDir, outputTreePath, Configs.numCores).run()
        treeutils.convertMafftGuideTree(outputTreePath, taxa)
    elif treeType.lower() == "clustal":
        Configs.log("Building Clustal Omega initial tree on {}..".format(context.sequencesPath))
        external_tools.runClustalOmegaGuideTree(context.sequencesPath, tempDir, outputTreePath, Configs.numCores).run()
    else:
        raise Exception("Guide tree {} not a file and not recognized..".format(treeType))

    time2 = time.time()
    Configs.log("Built initial tree on {} in {} sec..".format(context.sequencesPath, time2-time1))
    
    return outputTreePath

def buildInitialAlignment(sequences, tempDir, skeletonSize, initialAlignSize, outputAlignPath):
    skeletonPath = os.path.join(tempDir, "skeleton_sequences.txt")
    queriesPath = os.path.join(tempDir, "queries.txt") 
    hmmDir = os.path.join(tempDir, "skeleton_hmm")
    hmmPath = os.path.join(hmmDir, "hmm_model.txt")
    initialInsertPath = os.path.join(tempDir, "initial_insert_align.txt")
    if not os.path.exists(hmmDir):
        os.makedirs(hmmDir)
    
    if initialAlignSize is None or initialAlignSize > len(sequences):
        initialAlignSize = len(sequences)
    
    skeletonTaxa, remainingTaxa = decomposer.chooseSkeletonTaxa(sequences, skeletonSize)
    additional = initialAlignSize-skeletonSize
    random.shuffle(remainingTaxa)
    remainingTaxa, unusedTaxa = remainingTaxa[:additional], remainingTaxa[additional:]
    
    sequenceutils.writeFasta(sequences, skeletonPath, skeletonTaxa)
    external_tools.runMafft(skeletonPath, None, tempDir, outputAlignPath, Configs.numCores).run()
    
    if len(remainingTaxa) > 0:
        sequenceutils.writeFasta(sequences, queriesPath, remainingTaxa)    
        hmmutils.buildHmmOverAlignment(outputAlignPath, hmmPath).run()
        hmmTasks = hmmutils.hmmAlignQueries(hmmPath, queriesPath)
        task.submitTasks(hmmTasks)
        for hmmTask in task.asCompleted(hmmTasks):
            hmmutils.mergeHmmAlignments([hmmTask.outputFile], outputAlignPath, includeInsertions=False)
            if Configs.graphBuildMethod == "initial":
                hmmutils.mergeHmmAlignments([hmmTask.outputFile], initialInsertPath, includeInsertions=True)
