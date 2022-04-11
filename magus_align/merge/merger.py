'''
Created on May 14, 2020

@author: Vlad
'''

import os
import time


from magus_align.merge.graph_build.graph_builder import buildGraph
from magus_align.merge.graph_cluster.clusterer import clusterGraph
from magus_align.merge.graph_trace.tracer import findTrace
from magus_align.merge.optimizer import optimizeTrace
from magus_align.merge.alignment_writer import writeAlignment
from magus_configuration import Configs


def mergeSubalignments(context):
    Configs.log("Merging {} subaligments..".format(len(context.subalignmentPaths)))
    time1 = time.time()  
    
    buildGraph(context)
    clusterGraph(context.graph)
    findTrace(context.graph)
    optimizeTrace(context.graph)    
    writeAlignment(context)
    
    time2 = time.time()  
    Configs.log("Merged {} subalignments into {} in {} sec..".format(len(context.subalignmentPaths), context.outputFile, time2-time1))
