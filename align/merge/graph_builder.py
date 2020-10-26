'''
Created on Apr 14, 2020

@author: Vlad
'''

import os
import time
import random
import threading

from align.merge.alignment_graph import AlignmentGraph
from helpers import sequenceutils, hmmutils, tasks
from configuration import Configs
from tools import external_tools

def buildGraph(task):
    time1 = time.time() 
    
    task.graph = AlignmentGraph(task)
    
    if Configs.graphBuildMethod == "mafft":
        generateMafftBackbones(task)
        
    task.awaitSubalignments()
    task.graph.loadSubalignments(task.subalignmentPaths)
    
    if os.path.exists(task.graph.graphPath):
        task.graph.readGraphFromFile(task.graph.graphPath)
    else:
        buildMatrix(task)
        task.graph.writeGraphToFile(task.graph.graphPath)
        
    time2 = time.time()
    Configs.log("Built the alignment graph in {} sec..".format(time2-time1))


def buildMatrix(task):
    task.graph.matrix = [{} for i in range(task.graph.matrixSize)]
    task.graph.matrixLock = threading.Lock()
    
    if task.backbonePaths is not None and len(task.backbonePaths) > 0:
        Configs.log("Using {} user-defined backbone files..".format(len(task.backbonePaths)))
        for backboneAlignment in task.backbonePaths:
            addAlignmentFileToGraph(task.graph, backboneAlignment)
            
    elif Configs.graphBuildMethod == "mafft":
        applyMafftBackbones(task)
        
    elif Configs.graphBuildMethod == "hmm":
        applyHmmBackbones(task.graph)
        
    elif Configs.graphBuildMethod == "initial":
        Configs.log("Using the initial decomposition alignment as the single backbone..")
        mafftAlignPath = os.path.join(task.workingDir, "decomposition", "initial_tree", "skeleton_align.txt")
        hmmAlignPath = os.path.join(task.workingDir, "decomposition", "initial_tree", "queries_align.txt")
        addAlignmentFileToGraph(task.graph, mafftAlignPath, hmmAlignPath)
    

def generateMafftBackbones(task):
    Configs.log("Requesting {} MAFFT backbones..".format(Configs.mafftRuns)) 
    numSubsets = len(task.subsetPaths)
    subsets  = [sequenceutils.readFromFasta(filePath, removeDashes=True) for filePath in  task.subsetPaths]
    minSubalignmentLength = min([len(subset) for subset in subsets])
    backboneSubsetSize = max(1, min(minSubalignmentLength, int(Configs.mafftSize/numSubsets)))
    
    for n in range(Configs.mafftRuns):
        unalignedFile = os.path.join(task.graph.workingDir, "backbone_{}_unalign.txt".format(n+1))
        alignedFile = os.path.join(task.graph.workingDir, "backbone_{}_mafft.txt".format(n+1))
        if not os.path.exists(unalignedFile):
            assignBackboneTaxa(subsets, backboneSubsetSize, unalignedFile)
        backboneTask = external_tools.buildMafftAlignment(unalignedFile, alignedFile)
        task.backboneTasks.append(backboneTask)
    tasks.submitTasks(task.backboneTasks)

def applyMafftBackbones(task):    
    if Configs.graphBuildHmmExtend:
        taskExtensions = {}
        for backboneTask in tasks.asCompleted(task.backboneTasks):
            hmmTasks = createHmmExtensionTasks(task.graph, backboneTask.outputFile)
            tasks.submitTasks(hmmTasks)
            taskExtensions[backboneTask] = hmmTasks
        
    for backboneTask in tasks.asCompleted(task.backboneTasks):
        if Configs.graphBuildHmmExtend:
            baseName = os.path.basename(backboneTask.outputFile)
            extensionAlignedFile = os.path.join(task.graph.workingDir, "hmmextension_{}".format(baseName))
            
            for t in tasks.asCompleted(taskExtensions[backboneTask]):
                hmmutils.mergeHmmAlignments([t.outputFile], extensionAlignedFile, includeInsertions=True)
            #tasks.awaitTasks(taskExtensions[backboneTask])
            #alignFiles = [t.outputFile for t in taskExtensions[backboneTask]]
            #hmmutils.combineHmmAlignments(alignFiles, extensionAlignedFile, includeInsertions = True)
        else:
            extensionAlignedFile = None
        addAlignmentFileToGraph(task.graph, backboneTask.outputFile, extensionAlignedFile)
    
def assignBackboneTaxa(subsets, backboneSubsetSize, unalignedFile):
    backbone = {}
    for subset in subsets:
        taxa = list(subset.keys())
        random.shuffle(taxa)
        for taxon in taxa[:backboneSubsetSize]:
            backbone[taxon] = subset[taxon]
    sequenceutils.writeFasta(backbone, unalignedFile) 

def createHmmExtensionTasks(graph, alignedFile):
    baseName = os.path.basename(alignedFile)
    hmmDir = os.path.join(graph.workingDir, "hmm_{}".format(baseName))
    extensionUnalignedFile = os.path.join(hmmDir, "queries.txt")
    hmmPath = os.path.join(hmmDir, "hmm_model.txt")
    if not os.path.exists(hmmDir):
        os.makedirs(hmmDir)
    
    backbone = sequenceutils.readFromFasta(alignedFile)        
    backboneExtension = {}
    for taxon in graph.unalignment:
        if not taxon in backbone:
            backboneExtension[taxon] = graph.unalignment[taxon]
                
    sequenceutils.writeFasta(backboneExtension, extensionUnalignedFile) 
    buildTask = hmmutils.buildHmmOverAlignment(alignedFile, hmmPath)
    buildTask.run()
    alignTasks = hmmutils.hmmAlignQueries(hmmPath, extensionUnalignedFile)
    return alignTasks
    
def addAlignmentFileToGraph(graph, alignedFile, extensionAlignedFile = None):
    Configs.log("Feeding backbone {} to the graph..".format(alignedFile))
    backboneAlign = sequenceutils.readFromFasta(alignedFile)  
    alignmentLength = len(next(iter(backboneAlign.values())).seq)   
    
    if extensionAlignedFile is not None:
        extensionAlign = sequenceutils.readFromFasta(extensionAlignedFile)
        backboneAlign.update(extensionAlign)
        Configs.log("Applied HMM-extension to the backbone alignment from {}".format(alignedFile))
    
    alignmap = backboneToAlignMap(graph, backboneAlign, alignmentLength)
    Configs.log("Constructed backbone alignment map from {}".format(alignedFile))
    
    with graph.matrixLock:
        for l in range(alignmentLength):  
            for a, avalue in alignmap[l].items():
                for b, bvalue in alignmap[l].items():
                    
                    if Configs.graphBuildRestrict:
                        asub, apos = graph.matSubPosMap[a] 
                        bsub, bpos = graph.matSubPosMap[b] 
                        if asub == bsub and apos != bpos:
                            continue
                                
                    graph.matrix[a][b] = graph.matrix[a].get(b,0) + avalue * bvalue         
    Configs.log("Fed backbone {} to the graph.".format(alignedFile))

def backboneToAlignMap(graph, backboneAlign, alignmentLength):
    alignmap = [{} for i in range(alignmentLength)]
    t = 0
    
    for taxon in backboneAlign:
        subsetIdx = graph.taxonSubalignmentMap[taxon]
        subsetseq = graph.subalignments[subsetIdx][taxon].seq
        unalignedseq = graph.unalignment[taxon].seq     
        backboneseq = backboneAlign[taxon].seq   
        
        i = 0
        posarray = [0] * len(unalignedseq)
        for n in range(len(subsetseq)):
            if subsetseq[n] == unalignedseq[i]:
                posarray[i] = n 
                i = i + 1
                if i == len(unalignedseq):
                    break
        
        i = 0
        n = 0
        for c in backboneseq:
            if i == len(unalignedseq):
                break
            if c == unalignedseq[i]:
                position = int(graph.subsetMatrixIdx[subsetIdx] + posarray[i])
                alignmap[n][position] = alignmap[n].get(position, 0) + 1
            if c.upper() == unalignedseq[i]:
                i = i + 1
            if c == c.upper() and c != '.':
                n = n + 1
                
        t = t + 1
            
    return alignmap

def applyHmmBackbones(graph):
    numSubsets = len(graph.subalignments)
    Configs.log("Requesting {} HMM backbones..".format(numSubsets)) 
    
    subsetTasks = []
    for n in range(numSubsets):
        hmmTasks = createHmmExtensionTasks(graph, graph.subsetPaths[n])
        tasks.submitTasks(hmmTasks)
        subsetTasks.append(hmmTasks)
        
    for n, hmmTasks in enumerate(subsetTasks):
        tasks.awaitTasks(hmmTasks)
        baseName = os.path.basename(graph.subsetPaths[n])
        extensionAlignedFile = os.path.join(graph.workingDir, "extension_{}".format(baseName))
        hmmutils.combineHmmAlignments(tasks, extensionAlignedFile, includeInsertions = True)
        addAlignmentFileToGraph(graph, graph.subsetPaths[n], extensionAlignedFile)     
    