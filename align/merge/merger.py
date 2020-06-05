'''
Created on May 14, 2020

@author: Vlad
'''

import os

from align.merge.library_graph import LibraryGraph
from align.merge.graph_build import buildGraph
from align.merge.graph_cluster import clusterGraph
from align.merge.graph_order import orderGraph

def mergeSubalignments(workingDir, subalignmentPaths, outputPath):
    baseName = os.path.splitext(os.path.basename(outputPath))[0]
    mergingDir = os.path.join(workingDir, "merging_{}".format(baseName))
    if not os.path.exists(mergingDir):
        os.makedirs(mergingDir)
    
    graph = LibraryGraph(mergingDir)
    graph.loadSubalignments(subalignmentPaths)    
    buildGraph(graph)
    clusterGraph(graph)
    orderGraph(graph)
    graph.clustersToAlignment(outputPath)