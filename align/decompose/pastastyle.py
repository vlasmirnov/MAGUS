'''
Created on May 29, 2020

@author: Vlad
'''

import os
import shutil
import time

from align.decompose import decomposer
from helpers import sequenceutils, hmmutils
from tools import external_tools
from configuration import Configs


def buildPastaInitialTree(subsetsDir, sequencesPath):
    tempDir = os.path.join(subsetsDir, "initial_tree")
    outputTreePath = os.path.join(tempDir, "initial_tree.tre")
    outputAlignPath = os.path.join(tempDir, "initial_align.txt")
    
    if os.path.exists(outputTreePath):
        return outputTreePath, outputAlignPath
    if os.path.exists(tempDir):
        shutil.rmtree(tempDir)
    os.makedirs(tempDir)
    
    Configs.log("Building PASTA-style initial tree on {} with skeleton size {}..".format(sequencesPath, Configs.decompositionSkeletonSize))
    time1 = time.time() 
    
    skeletonPath = os.path.join(tempDir, "skeleton_sequences.txt")
    skeletonAlignPath = os.path.join(tempDir, "skeleton_align.txt") 
    queriesPath = os.path.join(tempDir, "queries.txt") 
    hmmDir = os.path.join(tempDir, "skeleton_hmm")
    hmmAlignPath = os.path.join(tempDir, "queries_align.txt")
    
    sequences = sequenceutils.readFromFasta(sequencesPath, True)
    skeletonTaxa, remainingTaxa = decomposer.chooseSkeletonTaxa(sequences, Configs.decompositionSkeletonSize)
    sequenceutils.writeFasta(sequences, skeletonPath, skeletonTaxa)
    sequenceutils.writeFasta(sequences, queriesPath, remainingTaxa)
    
    external_tools.runMafft(skeletonPath, None, tempDir, skeletonAlignPath, Configs.numCores)
    hmmutils.buildHmmOverAlignment(hmmDir, skeletonAlignPath)
    hmmutils.hmmAlignQueries(hmmDir, queriesPath, hmmAlignPath)
    initialAlign = sequenceutils.readFromStockholm(hmmAlignPath)
    skeletonAlign = sequenceutils.readFromFasta(skeletonAlignPath)
    initialAlign.update(skeletonAlign)
  
    sequenceutils.writeFasta(initialAlign, outputAlignPath)    
    external_tools.runFastTree(outputAlignPath, tempDir, outputTreePath)
    
    time2 = time.time()
    Configs.log("Built initial tree on {} in {} sec..".format(sequencesPath, time2-time1))
    
    return outputTreePath, outputAlignPath

    