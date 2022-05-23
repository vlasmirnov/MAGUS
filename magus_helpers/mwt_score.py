'''
Created on May 23, 2022

@author: Vlad
'''

import os
import argparse
from magus_configuration import Configs
from magus_align.alignment_context import AlignmentContext
from magus_align.merge.alignment_graph import AlignmentGraph
from magus_helpers import sequenceutils
from fileinput import filename


def computeClusteringScore(graph, clusters, cost = False):
        cutCost = 0
        nodeClusters = {}
        
        for n, cluster in enumerate(clusters):
            for a in cluster:
                nodeClusters[a] = n
                
        clusterCounter = len(clusters)
        for a in range(graph.matrixSize):
            if a not in nodeClusters:
                nodeClusters[a] = clusterCounter
                clusterCounter = clusterCounter + 1
                 
        for a in range(graph.matrixSize):
            asub, apos = graph.matSubPosMap[a]
            for b, value in graph.matrix[a].items():
                bsub, bpos = graph.matSubPosMap[b] 
                if asub != bsub and (nodeClusters[a] != nodeClusters[b]) == cost:
                    cutCost = cutCost + value
    
        return int(cutCost/2) 
    
def computeMWTScoreFromFiles(subalignmentPaths, matrixPath, outputPath, mode):
    #context = AlignmentContext(workingDir = os.path.dirname(outputPath), subalignmentPaths = subalignmentPaths)
    context = AlignmentContext(workingDir = os.path.dirname(os.path.dirname(matrixPath)), subalignmentPaths = subalignmentPaths)
    context.subsetPaths = context.subalignmentPaths
    context.initializeSequences()
    context.graph = AlignmentGraph(context)
    context.graph.graphPath = matrixPath
    context.graph.initializeMatrix()
    context.graph.readGraphFromFile(context.graph.graphPath)
    getClustersFromOutputFile(context, outputPath)
    print("MWT {}: {}".format(mode, computeClusteringScore(context.graph, context.graph.clusters, mode.lower() == "cost")))
    
def getClustersFromOutputFile(context, outputPath):
    outputAlign = sequenceutils.readFromFasta(outputPath, removeDashes = False)
    totalLen = len(next(iter(outputAlign.values())).seq)
    clusters = [[] for i in range(totalLen)]
    curPos = 0
    for i, path in enumerate(context.subalignmentPaths):
        subAlign = sequenceutils.readFromFasta(path, removeDashes = False)
        #curPos = 0
        for j in range(totalLen):
            for taxon in subAlign:
                if outputAlign[taxon].seq[j] not in ('-', '_'):
                    #clusters[j].append(graph.subsetMatrixIdx[i] + curPos)
                    clusters[j].append(curPos)
                    curPos = curPos + 1
                    break
    context.graph.clusters = clusters
    #multclusters = [c for c in clusters if len(c) > 1]
    #context.graph.clusters = multclusters
    #print(len(multclusters))
    print("Found {} clusters, {} / {} subalignment columns..".format(len(clusters), curPos, context.graph.matrixSize))
        
def main():   
    args = parseArgs()
    
    subalignmentPaths = []
    for p in args.subalignments:
        path = os.path.abspath(p)
        if os.path.isdir(path):
            for i in range(len(list(os.listdir(path)))):    
            #for filename in os.listdir(path):
                filename = os.path.join(path, "subalignment_subset_{}.txt".format(i+1))
                subalignmentPaths.append(os.path.join(path, filename))
        else:
            subalignmentPaths.append(path)
    matrixPath = os.path.abspath(args.library)
    outputPath = os.path.abspath(args.output)
    mode = args.mode
    computeMWTScoreFromFiles(subalignmentPaths, matrixPath, outputPath, mode)
    
def parseArgs():
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--subalignments", type=str, nargs="+",
                        help="Paths to input subalignment files", required=True, default=[])
    
    parser.add_argument("-l", "--library", type=str, 
                        help="Path to library graph file", required=True)
    
    parser.add_argument("-o", "--output", type=str,
                        help="Output alignment path", required=True)
    
    parser.add_argument("-m", "--mode", type=str,
                        help="cost or score", required=False, default = "cost")
    
       
    return parser.parse_args()

if __name__ == '__main__':
    main()