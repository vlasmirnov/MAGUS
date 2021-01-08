'''
Created on May 29, 2020

@author: Vlad
'''

import os
import shutil
import time
import random

from align.decompose import decomposer
from helpers import sequenceutils, hmmutils, treeutils
from tasks import task
from tools import external_tools
from configuration import Configs


def buildPastaInitialTree(context, subsetsDir):
    tempDir = os.path.join(subsetsDir, "initial_tree")
    outputTreePath = os.path.join(tempDir, "initial_tree.tre")
    outputAlignPath = os.path.join(tempDir, "initial_align.txt")
    
    if os.path.exists(outputTreePath):
        return outputTreePath, outputAlignPath
    if os.path.exists(tempDir):
        shutil.rmtree(tempDir)
    os.makedirs(tempDir)
    
    Configs.log("Building PASTA-style initial tree on {} with skeleton size {}..".format(context.sequencesPath, Configs.decompositionSkeletonSize))
    time1 = time.time() 

    taxa = list(context.unalignedSequences.keys())
    if len(taxa) >= 100000:
        external_tools.runMafftGuideTree(context.sequencesPath, tempDir, outputTreePath, Configs.numCores).run()
        treeutils.convertMafftGuideTree(outputTreePath, taxa)
    else:
        buildInitialAlignment(context.unalignedSequences, tempDir, Configs.decompositionSkeletonSize, None, outputAlignPath)
        external_tools.runFastTree(outputAlignPath, tempDir, outputTreePath, "fast").run()
    
    time2 = time.time()
    Configs.log("Built initial tree on {} in {} sec..".format(context.sequencesPath, time2-time1))
    
    return outputTreePath, outputAlignPath

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
