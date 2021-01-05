'''
Created on Apr 14, 2020

@author: Vlad
'''

import os
import time
import random
import gc

from align.merge.alignment_graph import AlignmentGraph
from helpers import sequenceutils, hmmutils
from tasks import task
from configuration import Configs
from tools import external_tools

def buildGraph(context):
    time1 = time.time() 
    
    context.graph = AlignmentGraph(context)
    context.initializeSequences()
    
    if os.path.exists(context.graph.graphPath):
        Configs.log("Found existing graph file {}".format(context.graph.graphPath))
    else:
        requestBackboneTasks(context)
    
    context.awaitSubalignments()
    context.graph.initializeMatrix()
    
    if os.path.exists(context.graph.graphPath):
        context.graph.readGraphFromFile(context.graph.graphPath)
    else:
        context.graph.initializeBackboneSequenceMapping()
        buildMatrix(context)
        context.graph.writeGraphToFile(context.graph.graphPath)
       
    gc.collect()   
    time2 = time.time()
    Configs.log("Built the alignment graph in {} sec..".format(time2-time1))

def requestBackboneTasks(context):
    if len(context.backbonePaths) > 0:
        Configs.log("Using {} user-defined backbone files..".format(len(context.backbonePaths)))
        context.graph.backbonePaths = context.backbonePaths
        for path in context.backbonePaths:
            context.graph.backboneTaxa.update(sequenceutils.readFromFasta(path))
    
    elif Configs.graphBuildMethod == "mafft":
        Configs.log("Using {} MAFFT backbones..".format(Configs.mafftRuns)) 
        requestMafftBackbones(context)
        
    elif Configs.graphBuildMethod == "subsethmm":
        Configs.log("Using {} HMM-extended subalignments as backbone files..".format(len(context.subalignmentPaths)))
        context.graph.backbonePaths = context.subalignmentPaths
        context.graph.backboneExtend.update(context.graph.backbonePaths)
        
    elif Configs.graphBuildMethod == "initial":
        Configs.log("Using the initial decomposition alignment as the single backbone..")
        initialAlignPath = os.path.join(context.workingDir, "decomposition", "initial_tree", "initial_align.txt")
        context.graph.backbonePaths = [initialAlignPath]

def requestMafftBackbones(context):
    numTaxa = max(1, int(Configs.mafftSize/len(context.subsetPaths)))

    for n in range(Configs.mafftRuns):
        unalignedFile = os.path.join(context.graph.workingDir, "backbone_{}_unalign.txt".format(n+1))
        alignedFile = os.path.join(context.graph.workingDir, "backbone_{}_mafft.txt".format(n+1))
        if os.path.exists(alignedFile):
            Configs.log("Existing backbone file found: {}".format(alignedFile))
            backbone = sequenceutils.readFromFasta(alignedFile)
            context.graph.backbonePaths.append(alignedFile)
        else:
            backbone = assignBackboneTaxa(context, numTaxa, unalignedFile)
            backboneTask = external_tools.buildMafftAlignment(unalignedFile, alignedFile)
            #backboneTask.submitTask()
            context.graph.backboneTasks.append(backboneTask)
            
        if Configs.graphBuildHmmExtend:
            context.graph.backboneExtend.add(alignedFile)
        else:
            context.graph.backboneTaxa.update(backbone)
    task.submitTasks(context.graph.backboneTasks)    
                
def buildMatrix(context):
    for backboneFile in context.graph.backbonePaths:
        addAlignmentFileToGraph(context.graph, backboneFile)
    
    for backboneTask in task.asCompleted(context.graph.backboneTasks):
        addAlignmentFileToGraph(context.graph, backboneTask.outputFile)
    
def assignBackboneTaxa(context, numTaxa, unalignedFile):
    backbone = {}
    for i, taxa in context.subsetTaxonMap.items():
        random.shuffle(taxa)
        for taxon in taxa[:numTaxa]:
            backbone[taxon] = context.unalignedSequences[taxon]
    sequenceutils.writeFasta(backbone, unalignedFile)
    return backbone
    
def addAlignmentFileToGraph(graph, alignedFile):
    Configs.log("Feeding backbone {} to the graph..".format(alignedFile))
    backboneAlign = sequenceutils.readFromFasta(alignedFile)  
    alignmentLength = len(next(iter(backboneAlign.values())).seq)   
    
    if alignedFile in graph.backboneExtend:
        extensionTasks = requestHmmExtensionTasks(graph, backboneAlign, alignedFile)
        task.submitTasks(extensionTasks)
        for extensionTask in task.asCompleted(extensionTasks):
            backboneAlign.update(sequenceutils.readFromStockholm(extensionTask.outputFile, includeInsertions=True))
    
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
        subsetIdx = graph.context.taxonSubsetMap[taxon]
        subsetseq = graph.backboneSubsetAlignment[taxon].seq
        unalignedseq = graph.context.unalignedSequences[taxon].seq     
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
        
def requestHmmExtensionTasks(graph, backbone, alignedFile):
    baseName = os.path.basename(alignedFile)
    hmmDir = os.path.join(graph.workingDir, "hmm_{}".format(baseName))
    extensionUnalignedFile = os.path.join(hmmDir, "queries.txt")
    hmmPath = os.path.join(hmmDir, "hmm_model.txt")
    if not os.path.exists(hmmDir):
        os.makedirs(hmmDir)
    
    backboneExtension = {}
    for taxon in graph.context.unalignedSequences:
        if not taxon in backbone:
            backboneExtension[taxon] = graph.context.unalignedSequences[taxon]
                
    sequenceutils.writeFasta(backboneExtension, extensionUnalignedFile) 
    buildTask = hmmutils.buildHmmOverAlignment(alignedFile, hmmPath)
    buildTask.run()
    alignTasks = hmmutils.hmmAlignQueries(hmmPath, extensionUnalignedFile)
    return alignTasks 
    