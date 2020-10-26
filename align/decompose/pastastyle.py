'''
Created on May 29, 2020

@author: Vlad
'''

import os
import shutil
import time
import random

from align.decompose import decomposer
from helpers import sequenceutils, hmmutils, treeutils, tasks
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
    
    align = sequenceutils.readFromFasta(sequencesPath)
    taxa = list(align.keys())
    if len(taxa) >= 100000:
        external_tools.runMafftGuideTree(sequencesPath, tempDir, outputTreePath, Configs.numCores).run()
        treeutils.convertMafftGuideTree(outputTreePath, taxa)
        
        #outputTreeMafftPath = os.path.join(tempDir, "initial_tree_mafft.tre")
        #external_tools.runMafftGuideTree(sequencesPath, tempDir, outputTreeMafftPath, Configs.numCores).run()
        #treeutils.convertMafftGuideTree(outputTreeMafftPath, taxa)
        
        #initialAlign, unused = buildInitialAlignment(sequencesPath, tempDir, Configs.decompositionSkeletonSize, None)
        #sequenceutils.writeFasta(initialAlign, outputAlignPath)    
        #external_tools.runFastTree(outputAlignPath, tempDir, outputTreePath, "fastest", intree=outputTreeMafftPath).run()
    else:
        initialAlign, unused = buildInitialAlignment(sequencesPath, tempDir, Configs.decompositionSkeletonSize, None)
        sequenceutils.writeFasta(initialAlign, outputAlignPath)    
        if len(initialAlign) >= 50000:
            external_tools.runFastTree(outputAlignPath, tempDir, outputTreePath, "fastest").run()
        elif len(initialAlign) >= 20000:
            external_tools.runFastTree(outputAlignPath, tempDir, outputTreePath, "faster").run()
        elif len(initialAlign) >= 5000:
            external_tools.runFastTree(outputAlignPath, tempDir, outputTreePath, "fast").run()
        else:
            external_tools.runFastTree(outputAlignPath, tempDir, outputTreePath, "normal").run()
    
    '''
    align = sequenceutils.readFromFasta(sequencesPath)
    taxa = list(align.keys())
    if len(taxa) >= 10000:
        external_tools.runMafftGuideTree(sequencesPath, tempDir, outputTreePath, Configs.numCores).run()
        treeutils.convertMafftGuideTree(outputTreePath, taxa)
    else:
        initialAlign, unused = buildInitialAlignment(sequencesPath, tempDir, Configs.decompositionSkeletonSize, None)
        sequenceutils.writeFasta(initialAlign, outputAlignPath)    
        if len(initialAlign) >= 10000:
            external_tools.runFastTree(outputAlignPath, tempDir, outputTreePath, "fastest").run()
        elif len(initialAlign) >= 5000:
            external_tools.runFastTree(outputAlignPath, tempDir, outputTreePath, "fast").run()
        else:
            external_tools.runFastTree(outputAlignPath, tempDir, outputTreePath, "normal").run()
    '''
    
    
    time2 = time.time()
    Configs.log("Built initial tree on {} in {} sec..".format(sequencesPath, time2-time1))
    
    return outputTreePath, outputAlignPath

def buildInitialAlignment(sequencesPath, tempDir, skeletonSize, initialAlignSize):
    skeletonPath = os.path.join(tempDir, "skeleton_sequences.txt")
    skeletonAlignPath = os.path.join(tempDir, "skeleton_align.txt") 
    queriesPath = os.path.join(tempDir, "queries.txt") 
    hmmDir = os.path.join(tempDir, "skeleton_hmm")
    hmmPath = os.path.join(hmmDir, "hmm_model.txt")
    hmmAlignPath = os.path.join(tempDir, "queries_align.txt")
    if not os.path.exists(hmmDir):
        os.makedirs(hmmDir)
    
    sequences = sequenceutils.readFromFasta(sequencesPath, True)
    if initialAlignSize is None or initialAlignSize > len(sequences):
        initialAlignSize = len(sequences)
    
    skeletonTaxa, remainingTaxa = decomposer.chooseSkeletonTaxa(sequences, skeletonSize)
    additional = initialAlignSize-skeletonSize
    random.shuffle(remainingTaxa)
    remainingTaxa, unusedTaxa = remainingTaxa[:additional], remainingTaxa[additional:]
    
    sequenceutils.writeFasta(sequences, skeletonPath, skeletonTaxa)
    external_tools.runMafft(skeletonPath, None, tempDir, skeletonAlignPath, Configs.numCores).run()
    initialAlign = sequenceutils.readFromFasta(skeletonAlignPath)
    
    if len(remainingTaxa) > 0:
        sequenceutils.writeFasta(sequences, queriesPath, remainingTaxa)    
        hmmutils.buildHmmOverAlignment(skeletonAlignPath, hmmPath).run()
        hmmTasks = hmmutils.hmmAlignQueries(hmmPath, queriesPath)
        tasks.submitTasks(hmmTasks)
        for task in tasks.asCompleted(hmmTasks):
            hmmutils.mergeHmmAlignments([task.outputFile], hmmAlignPath, includeInsertions=False)
        #hmmutils.combineHmmAlignments([t.outputFile for t in hmmTasks], hmmAlignPath, includeInsertions=False)
        hmmAlign = sequenceutils.readFromFasta(hmmAlignPath)
        initialAlign.update(hmmAlign)
    return initialAlign, unusedTaxa
    