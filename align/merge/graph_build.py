'''
Created on Apr 14, 2020

@author: Vlad
'''

import os
import time
import random
import concurrent.futures
import threading

from helpers import sequenceutils
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
    Configs.log("Launching {} backbones with {} workers..".format(Configs.mafftRuns, Configs.numCores))
        
    numSubsets = len(graph.subalignments)
    minSubalignmentLength = min([len(subset) for subset in graph.subalignments])
    backboneSubsetSize = max(1, min(minSubalignmentLength, int(Configs.mafftSize/numSubsets)))
    
    with concurrent.futures.ThreadPoolExecutor(max_workers = Configs.numCores) as executor:
        jobs = {executor.submit(addMafftBackboneToGraph, graph, backboneSubsetSize, n+1) : n+1
                   for n in range(Configs.mafftRuns)}
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
    if not os.path.exists(unalignedFile):
        backbone = {}
        for subset in graph.subalignments:
            taxa = list(subset.keys())
            random.shuffle(taxa)
            
            for taxon in taxa[:backboneSubsetSize]:
                backbone[taxon] = subset[taxon]

        sequenceutils.writeFasta(backbone, unalignedFile) 
    
    external_tools.buildMafftAlignment(unalignedFile, alignedFile)
    
    Configs.log("Feeding backbone {} to the graph..".format(backboneIndex))
    addAlignmentFileToGraph(alignedFile, graph)                            
    Configs.log("Fed backbone {} to the graph.".format(backboneIndex))

def addAlignmentFileToGraph(alignedFile, graph):
    backboneAlign = sequenceutils.readFromFasta(alignedFile)        
    alignvector = backboneToAlignVector(graph, backboneAlign, graph.subsetMatrixIdx)
    lenAlignment = len(next(iter(backboneAlign.values())).seq) 
    
    with graph.matrixLock:
        for l in range(lenAlignment):  
            for k1 in range(len(backboneAlign)):
                if alignvector[k1][l] != -1:
                    for i1 in range(len(backboneAlign)):
                        if alignvector[i1][l] != -1:
                            a, b = alignvector[k1][l], alignvector[i1][l]
                            graph.matrix[a][b] = graph.matrix[a].get(b,0) + 1
             

def backboneToAlignVector(graph, backboneAlign, subsetMatrixIdx):
    lenAlignment = len(next(iter(backboneAlign.values())).seq)    
    #alignvector = numpy.full((len(joinsetAlign), lenAlignment), -1)
    alignvector = [[-1] * lenAlignment for i in range(len(backboneAlign))]
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
        for n in range(len(backboneseq)):
            if i < len(unalignedseq) and backboneseq[n] == unalignedseq[i]:
                alignvector[t][n] = int(subsetMatrixIdx[subsetIdx] + posarray[i])
                i = i + 1
        t = t + 1
            
    return alignvector

        
    