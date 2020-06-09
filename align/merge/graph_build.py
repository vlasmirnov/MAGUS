'''
Created on Apr 14, 2020

@author: Vlad
'''

import os
import time
import random
import concurrent.futures
import threading
import string

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
    Configs.log("Compiled backbone alignments and graph in {} sec..".format(time2-time1))
        
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
        
        if Configs.libraryGraphStrategy == "mafftbackbones":
            Configs.log("Launching {} MAFFT backbones with {} workers..".format(Configs.mafftRuns, Configs.numCores)) 
            minSubalignmentLength = min([len(subset) for subset in graph.subalignments])
            backboneSubsetSize = max(1, min(minSubalignmentLength, int(Configs.mafftSize/numSubsets)))
            jobs = {executor.submit(addMafftBackboneToGraph, graph, backboneSubsetSize, n+1) : n+1
                       for n in range(Configs.mafftRuns)}
            
        elif Configs.libraryGraphStrategy == "hmmbackbones":
            Configs.log("Launching {} HMM backbones with {} workers..".format(numSubsets, Configs.numCores)) 
            jobs = {executor.submit(addHmmBackboneToGraph, graph, n+1) : n+1
                       for n in range(numSubsets)}
        
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
        
        if Configs.libraryGraphMafftMerge:
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
    
    
    Configs.log("Feeding backbone {} to the graph..".format(backboneIndex))
    addAlignmentFileToGraph(alignedFile, graph)                            
    Configs.log("Fed backbone {} to the graph.".format(backboneIndex))

def addAlignmentFileToGraph(alignedFile, graph):
    backboneAlign = sequenceutils.readFromFasta(alignedFile)  
    alignmentLength = len(next(iter(backboneAlign.values())).seq)         
    alignvector = backboneToAlignVector(graph, backboneAlign, alignmentLength)
    
    with graph.matrixLock:
        for l in range(alignmentLength):  
            for k1 in range(len(backboneAlign)):
                if alignvector[k1][l] != -1:
                    for i1 in range(len(backboneAlign)):
                        if alignvector[i1][l] != -1:                            
                            a, b = alignvector[k1][l], alignvector[i1][l]
                            
                            if Configs.libraryGraphRestrict:
                                asub, apos = graph.matSubPosMap[a] 
                                bsub, bpos = graph.matSubPosMap[b] 
                                if asub == bsub and apos != bpos:
                                    continue
                                
                            graph.matrix[a][b] = graph.matrix[a].get(b,0) + 1
             

def backboneToAlignVector(graph, backboneAlign, alignmentLength):
    alignvector = [[-1] * alignmentLength for i in range(len(backboneAlign))]
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
                alignvector[t][n] = int(graph.subsetMatrixIdx[subsetIdx] + posarray[i])       
            if c.upper() == unalignedseq[i]:
                i = i + 1
            if c == c.upper() and c != '.':
                n = n + 1
                
        t = t + 1
            
    return alignvector

def addHmmBackboneToGraph(graph, subsetIdx):
    unalignedFile = os.path.join(graph.workingDir, "backbone_{}_unalign.txt".format(subsetIdx))
    hmmDir = os.path.join(graph.workingDir, "backbone_{}_hmm".format(subsetIdx))
    hmmAlignmentPath = os.path.join(hmmDir, "hmm_align_result.txt")
    
    if not os.path.exists(hmmAlignmentPath):
        backbone = {}
        for i, subset in enumerate(graph.subalignments):
            if i == subsetIdx-1:
                continue
            for taxon in subset:
                backbone[taxon] = graph.unalignment[taxon]
                    
        sequenceutils.writeFasta(backbone, unalignedFile) 
        hmmutils.buildHmmOverAlignment(hmmDir, graph.subsetPaths[subsetIdx-1])
        hmmutils.hmmAlignQueries(hmmDir, unalignedFile)
    
    Configs.log("Feeding backbone {} to the graph..".format(subsetIdx))
    addHmmAlignmentFileToGraphLinear(hmmAlignmentPath, subsetIdx-1, graph)                            
    Configs.log("Fed backbone {} to the graph.".format(subsetIdx))

def addHmmAlignmentFileToGraphLinear(alignedFile, subsetIdx, graph):
    subalignment = graph.subalignments[subsetIdx]
    hmmAlignment = sequenceutils.readFromStockholm(alignedFile, includeInsertions = True)
    alignmentLength = len(next(iter(subalignment.values())).seq) 
    alignvector = backboneToAlignVector(graph, hmmAlignment, alignmentLength)
        
    with graph.matrixLock:
        for l in range(alignmentLength):  
            b = int(graph.subsetMatrixIdx[subsetIdx] + l)
            graph.matrix[b][b] = graph.matrix[b].get(b,0) + int(2*len(subalignment))  
            for k1 in range(len(hmmAlignment)):
                if alignvector[k1][l] != -1:
                    a = alignvector[k1][l]                        
                    graph.matrix[a][b] = graph.matrix[a].get(b,0) + 1    
                    graph.matrix[b][a] = graph.matrix[b].get(a,0) + 1    
                        
   
  