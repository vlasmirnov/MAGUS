'''
Created on May 14, 2020

@author: Vlad
'''

import os

from align.merge.alignment_graph import AlignmentGraph
from align.merge.graph_builder import buildGraph
from align.merge.graph_cluster.clusterer import clusterGraph
from align.merge.graph_trace.tracer import findTrace
from align.merge.optimizer import optimizeTrace

def mergeSubalignments(workingDir, subalignmentPaths, outputPath):
    baseName = os.path.splitext(os.path.basename(outputPath))[0]
    mergingDir = os.path.join(workingDir, "merging_{}".format(baseName))
    if not os.path.exists(mergingDir):
        os.makedirs(mergingDir)
    
    graph = AlignmentGraph(mergingDir)
    graph.loadSubalignments(subalignmentPaths)    
    buildGraph(graph)
    clusterGraph(graph)
    findTrace(graph)
    optimizeTrace(graph)
    graph.clustersToAlignment(outputPath)
    
    
