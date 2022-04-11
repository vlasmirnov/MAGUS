'''
Created on Apr 14, 2020

@author: Vlad
'''

import os
import time
import random

from magus_align.merge.alignment_graph import AlignmentGraph
from magus_helpers import sequenceutils, hmmutils
from magus_tasks import task
from magus_configuration import Configs
from magus_tools import external_tools

'''
Building a MAGUS alignment graph from backbone alignments. 
The graph is a sparse matrix, stored as a weighted adjacency list.
Backbone alignment tasks are run in parallel, and can begin before the subalignments are finished.
When subalignments finish, we can initialize the alignment graph. 
Backbones are added to this graph as they complete.
Graph is written to file in MCL-compliant format.
'''

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
        context.initializeBackboneSequenceMapping()
        buildMatrix(context)
        context.graph.writeGraphToFile(context.graph.graphPath)
       
    time2 = time.time()
    Configs.log("Built the alignment graph in {} sec..".format(time2-time1))

def requestBackboneTasks(context):
    if len(context.backbonePaths) > 0:
        Configs.log("Using {} user-defined backbone files..".format(len(context.backbonePaths)))
        context.backbonePaths = context.backbonePaths
        for path in context.backbonePaths:
            context.backboneTaxa.update(sequenceutils.readFromFasta(path))
    
    elif Configs.graphBuildMethod == "mafft":
        Configs.log("Using {} MAFFT backbones..".format(Configs.mafftRuns)) 
        requestMafftBackbones(context)
        
    elif Configs.graphBuildMethod == "subsethmm":
        Configs.log("Using {} HMM-extended subalignments as backbone files..".format(len(context.subalignmentPaths)))
        context.backbonePaths = context.subalignmentPaths
        context.backboneExtend.update(context.backbonePaths)
        
    elif Configs.graphBuildMethod == "initial":
        Configs.log("Using the initial decomposition alignment as the single backbone..")
        initialAlignPath = os.path.join(context.workingDir, "decomposition", "initial_tree", "initial_insert_align.txt")
        context.backbonePaths = [initialAlignPath]
    
    if not Configs.constrain and Configs.graphBuildMethod != "subsethmm":
        context.backbonePaths.extend(context.subalignmentPaths)

def requestMafftBackbones(context):
    numTaxa = max(1, int(Configs.mafftSize/len(context.subsetPaths)))

    for n in range(Configs.mafftRuns):
        unalignedFile = os.path.join(context.graph.workingDir, "backbone_{}_unalign.txt".format(n+1))
        alignedFile = os.path.join(context.graph.workingDir, "backbone_{}_mafft.txt".format(n+1))
        if os.path.exists(alignedFile):
            Configs.log("Existing backbone file found: {}".format(alignedFile))
            backbone = sequenceutils.readFromFasta(alignedFile)
            context.backbonePaths.append(alignedFile)
        else:
            backbone = assignBackboneTaxa(context, numTaxa, unalignedFile)
            backboneTask = external_tools.buildMafftAlignment(unalignedFile, alignedFile)
            context.backboneTasks.append(backboneTask)
            
        if Configs.graphBuildHmmExtend:
            context.backboneExtend.add(alignedFile)
        else:
            context.backboneTaxa.update(backbone)
    task.submitTasks(context.backboneTasks)    
                
def buildMatrix(context):
    addedBackbones = set()
    for backboneTask in task.asCompleted(context.backboneTasks):
        addAlignmentFileToGraph(context, backboneTask.outputFile)
        addedBackbones.add(backboneTask.outputFile)
    
    for backboneFile in context.backbonePaths:
        if backboneFile not in addedBackbones:
            addAlignmentFileToGraph(context, backboneFile)  
    
def assignBackboneTaxa(context, numTaxa, unalignedFile):
    backbone = {}
    for subset in context.subsets:
        random.shuffle(subset)
        for taxon in subset[:numTaxa]:
            backbone[taxon] = context.unalignedSequences[taxon]
    sequenceutils.writeFasta(backbone, unalignedFile)
    return backbone
    
def addAlignmentFileToGraph(context, alignedFile):
    Configs.log("Feeding backbone {} to the graph..".format(alignedFile))
    backboneAlign = sequenceutils.readFromFasta(alignedFile)  
    alignmentLength = len(next(iter(backboneAlign.values())).seq)
     
    if alignedFile in context.backboneExtend:
        extensionTasks = requestHmmExtensionTasks(context, backboneAlign, alignedFile)
        task.submitTasks(extensionTasks)
        for extensionTask in task.asCompleted(extensionTasks):
            backboneAlign.update(sequenceutils.readFromStockholm(extensionTask.outputFile, includeInsertions=True))
    
    alignmap = backboneToAlignMap(context, backboneAlign, alignmentLength)
    Configs.log("Constructed backbone alignment map from {}".format(alignedFile))
    
    graph = context.graph  
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

def backboneToAlignMap(context, backboneAlign, alignmentLength):
    alignmap = [{} for i in range(alignmentLength)]
    t = 0
    
    for taxon in backboneAlign:
        subsetIdx = context.taxonSubalignmentMap[taxon]
        subsetseq = context.backboneSubalignment[taxon].seq
        unalignedseq = context.unalignedSequences[taxon].seq     
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
                position = int(context.graph.subsetMatrixIdx[subsetIdx] + posarray[i])
                alignmap[n][position] = alignmap[n].get(position, 0) + 1
            if c.upper() == unalignedseq[i]:
                i = i + 1
            if c == c.upper() and c != '.':
                n = n + 1
                
        t = t + 1
            
    return alignmap
        
def requestHmmExtensionTasks(context, backbone, alignedFile):
    baseName = os.path.basename(alignedFile)
    hmmDir = os.path.join(context.graph.workingDir, "hmm_{}".format(baseName))
    extensionUnalignedFile = os.path.join(hmmDir, "queries.txt")
    hmmPath = os.path.join(hmmDir, "hmm_model.txt")
    if not os.path.exists(hmmDir):
        os.makedirs(hmmDir)
    
    backboneExtension = {}
    for taxon in context.unalignedSequences:
        if not taxon in backbone:
            backboneExtension[taxon] = context.unalignedSequences[taxon]
                
    sequenceutils.writeFasta(backboneExtension, extensionUnalignedFile) 
    buildTask = hmmutils.buildHmmOverAlignment(alignedFile, hmmPath)
    buildTask.run()
    alignTasks = hmmutils.hmmAlignQueries(hmmPath, extensionUnalignedFile)
    return alignTasks 
    