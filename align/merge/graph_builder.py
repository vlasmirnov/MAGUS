'''
Created on Apr 14, 2020

@author: Vlad
'''

import os
import time
import random
import concurrent.futures
import threading

from helpers import sequenceutils, hmmutils
from configuration import Configs
from tools import external_tools

def buildGraph(graph):
    time1 = time.time() 
    
    graph.graphPath = os.path.join(graph.workingDir, "library_graph.txt")
    
    if os.path.exists(graph.graphPath):
        graph.readGraphFromFile(graph.graphPath)
    else:
        buildMatrix(graph)
        graph.writeGraphToFile(graph.graphPath)
        
    time2 = time.time()
    Configs.log("Built the alignment graph in {} sec..".format(time2-time1))
        
def buildMatrix(graph):
    graph.matrix = [{} for i in range(graph.matrixSize)]
    graph.matrixLock = threading.Lock()
    
    if Configs.backbonePaths is not None and len(Configs.backbonePaths) > 0:
        for backboneAlignment in Configs.backbonePaths:
            Configs.log("Feeding backbone file {} to the graph..".format(backboneAlignment))
            addAlignmentFileToGraph(backboneAlignment, graph)
            Configs.log("Fed backbone file {} to the graph.".format(backboneAlignment))
            
    else:
        buildBackboneAlignmentsAndMatrix(graph)
        
def buildBackboneAlignmentsAndMatrix(graph):
    numSubsets = len(graph.subalignments)
    
    with concurrent.futures.ThreadPoolExecutor(max_workers = Configs.numCores) as executor:
        
        if Configs.graphBuildMethod in ["mafft", "mafftmerge"]:
            Configs.log("Launching {} MAFFT backbones with {} workers..".format(Configs.mafftRuns, Configs.numCores)) 
            minSubalignmentLength = min([len(subset) for subset in graph.subalignments])
            backboneSubsetSize = max(1, min(minSubalignmentLength, int(Configs.mafftSize/numSubsets)))
            jobs = {executor.submit(addMafftBackboneToGraph, graph, backboneSubsetSize, n+1) : n+1
                       for n in range(Configs.mafftRuns)}
            
        elif Configs.graphBuildMethod == "hmm":
            Configs.log("Launching {} HMM backbones with {} workers..".format(numSubsets, Configs.numCores)) 
            jobs = {executor.submit(addHmmBackboneToGraph, graph, n+1) : n+1
                       for n in range(numSubsets)}
        
        elif Configs.graphBuildMethod == "initial":
            Configs.log("Using the initial decomposition alignment as the single backbone..")
            subsetsDir = os.path.dirname(graph.subsetPaths[0])
            mafftAlignPath = os.path.join(subsetsDir, "initial_tree", "skeleton_align.txt")
            hmmAlignPath = os.path.join(subsetsDir, "initial_tree", "queries_align.txt")
            jobs = {executor.submit(addAlignmentFileToGraph, mafftAlignPath, graph, hmmAlignPath) : 1}
        
        for job in concurrent.futures.as_completed(jobs):
            try:
                job.result()
            except Exception as exc:
                Configs.log("Worker for backbone {} threw an exception:\n{}".format(jobs[job], exc))
                raise
                        
    Configs.log("All workers done..")

def addMafftBackboneToGraph(graph, backboneSubsetSize, backboneIndex):
    unalignedFile = os.path.join(graph.workingDir, "backbone_{}_unalign.txt".format(backboneIndex))
    alignedFile = os.path.join(graph.workingDir, "backbone_{}_mafft.txt".format(backboneIndex))
    extensionAlignedFile = None
    if not os.path.exists(alignedFile):
        backbone = {}
        backboneTaxa = []
        for subset in graph.subalignments:
            taxa = list(subset.keys())
            random.shuffle(taxa)
            
            for taxon in taxa[:backboneSubsetSize]:
                backbone[taxon] = subset[taxon]
                backboneTaxa.append(taxon)

        sequenceutils.writeFasta(backbone, unalignedFile, backboneTaxa) 
        
        if Configs.graphBuildMethod == "mafftmerge":
            subtableFile = os.path.join(graph.workingDir, "backbone_{}_subtable.txt".format(backboneIndex))
            with open(subtableFile, 'w') as textFile:
                curSubset = None  
                for i, taxon in enumerate(backboneTaxa):
                    taxonSubset = graph.taxonSubalignmentMap[taxon]
                    if curSubset is not None and taxonSubset != curSubset:
                        textFile.write('\n')
                    curSubset = taxonSubset
                    textFile.write(str(i+1) + ' ')

            external_tools.buildMafftAlignment(unalignedFile, alignedFile, subtableFile)
        else:
            external_tools.buildMafftAlignment(unalignedFile, alignedFile)
    
        if Configs.graphBuildHmmExtend:
            extensionUnalignedFile = os.path.join(graph.workingDir, "backbone_{}_extension_unalign.txt".format(backboneIndex))
            hmmDir = os.path.join(graph.workingDir, "backbone_{}_hmm".format(backboneIndex))
            extensionAlignedFile = os.path.join(graph.workingDir, "backbone_{}_extension_align.txt".format(backboneIndex))
            
            backboneExtension = {}
            for taxon in graph.unalignment:
                if not taxon in backbone:
                    backboneExtension[taxon] = graph.unalignment[taxon]
                        
            sequenceutils.writeFasta(backboneExtension, extensionUnalignedFile) 
            hmmutils.buildHmmOverAlignment(hmmDir, alignedFile)
            hmmutils.hmmAlignQueries(hmmDir, extensionUnalignedFile, extensionAlignedFile)
    
    Configs.log("Feeding backbone {} to the graph..".format(backboneIndex))
    addAlignmentFileToGraph(alignedFile, graph, extensionAlignedFile)                            
    Configs.log("Fed backbone {} to the graph.".format(backboneIndex))

def addAlignmentFileToGraph(alignedFile, graph, extensionAlignedFile = None):
    backboneAlign = sequenceutils.readFromFasta(alignedFile)  
    alignmentLength = len(next(iter(backboneAlign.values())).seq)   
    
    if extensionAlignedFile is not None:
        extensionAlign = sequenceutils.readFromStockholm(extensionAlignedFile, includeInsertions = True)
        backboneAlign.update(extensionAlign)
    
    alignmap = backboneToAlignMap(graph, backboneAlign, alignmentLength)
    
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

def addHmmBackboneToGraph(graph, subsetIdx):
    subalignmentFile = graph.subsetPaths[subsetIdx-1]
    unalignedFile = os.path.join(graph.workingDir, "backbone_{}_unalign.txt".format(subsetIdx))
    hmmDir = os.path.join(graph.workingDir, "backbone_{}_hmm".format(subsetIdx))
    hmmAlignmentPath = os.path.join(graph.workingDir, "backbone_{}_align.txt".format(subsetIdx))
    
    if not os.path.exists(hmmAlignmentPath):
        backbone = {}
        for i, subset in enumerate(graph.subalignments):
            if i == subsetIdx-1:
                continue
            for taxon in subset:
                backbone[taxon] = graph.unalignment[taxon]
                    
        sequenceutils.writeFasta(backbone, unalignedFile) 
        hmmutils.buildHmmOverAlignment(hmmDir, subalignmentFile)
        hmmutils.hmmAlignQueries(hmmDir, unalignedFile, hmmAlignmentPath)
    
    Configs.log("Feeding backbone {} to the graph..".format(subsetIdx))
    addAlignmentFileToGraph(subalignmentFile, graph, hmmAlignmentPath)                                 
    Configs.log("Fed backbone {} to the graph.".format(subsetIdx))
