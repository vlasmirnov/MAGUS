'''
Created on Apr 14, 2020

@author: Vlad
'''

import time
import argparse
import sys

from gcm.library_graph import LibraryGraph
from gcm.graph_build import buildGraph
from gcm.graph_cluster import clusterGraph
from gcm.graph_order import orderGraph
from configuration import buildConfigs, Configs


def main(args):   
    startTime = time.time()
    buildConfigs(args)
    Configs.log("GCM was run with: {}".format(" ".join(sys.argv)))
    
    graph = LibraryGraph()
    graph.loadSubalignments(Configs.subsetPaths)    
    buildGraph(graph)
    clusterGraph(graph)
    orderGraph(graph)
    graph.clustersToAlignment(Configs.outputPath)
    
    endTime = time.time()
    Configs.log("GCM finished in {} seconds..".format(endTime-startTime))
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str,
                        help="Path to working directory", required=False, default=None)
    
    #parser.add_argument("-i", "--sequences", type=str,
    #                    help="Path to input unaligned sequences", required=False, default=None)

    parser.add_argument("-s", "--subalignments", type=str, nargs="+",
                        help="Paths to input subalignment files", required=True)

    parser.add_argument("-o", "--output", type=str,
                        help="Output alignment path",
                        required=True)
    
    #parser.add_argument("-t", "--guidetree", type=str,
    #                    help="guide tree for merge weights",
    #                    required=False, default=None)
    
    parser.add_argument("-r", "--mafftruns", type=int,
                        help="Number of MAFFT runs", required=False, default=10)
    
    parser.add_argument("-m", "--mafftsize", type=int,
                        help="Maximum size of MAFFT alignments", required=False, default=200)

    main(parser.parse_args())