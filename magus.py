'''
Created on Apr 14, 2020

@author: Vlad
'''

import time
import argparse
import sys
import traceback

from magus_align.aligner import mainAlignmentTask
from magus_configuration import buildConfigs, Configs
from magus_tasks import manager

def main():   
    '''
    Resolve the args/configs, spin up the task manager (which deals with worker threads and handles parallelism), 
    and get started on the main alignment task. 
    '''
    
    startTime = time.time()
    args = parseArgs()
    buildConfigs(args)    
    Configs.log("MAGUS was run with: {}".format(" ".join(sys.argv)))
    
    try:
        manager.startTaskManager()
        mainAlignmentTask()
    except:
        Configs.error("MAGUS aborted with an exception..")
        Configs.error(traceback.format_exc())
    finally:
        manager.stopTaskManager()
    
    endTime = time.time()
    Configs.log("MAGUS finished in {} seconds..".format(endTime-startTime))
    
def parseArgs():
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str,
                        help="Path to working directory", required=False, default=None)
    
    parser.add_argument("-i", "--sequences", type=str,
                        help="Path to input unaligned sequences", required=False, default=None)
    
    parser.add_argument("-s", "--subalignments", type=str, nargs="+",
                        help="Paths to input subalignment files", required=False, default=[])
    
    parser.add_argument("-b", "--backbones", type=str, nargs="+",
                        help="Paths to input backbone alignment files", required=False, default=[])

    parser.add_argument("-o", "--output", type=str,
                        help="Output alignment path", required=True)
    
    parser.add_argument("-t", "--guidetree", type=str,
                        help="Guide tree for subset decomposition. fasttree (default), fasttree-noml, clustal, parttree, or path to user guide tree",
                        required=False, default="fasttree")

    parser.add_argument("-np", "--numprocs", type=int,
                        help="Number of processors to use (default: # cpus available)",
                        required=False, default=-1)
    
    parser.add_argument("--maxsubsetsize", type=int,
                        help="Maximum subset size for divide-and-conquer",
                        required=False, default=50)
    
    parser.add_argument("--maxnumsubsets", type=int,
                        help="Maximum number of subsets for divide-and-conquer",
                        required=False, default=25)
    
    parser.add_argument("--decompstrategy", type=str,
                        help="Initial decomposition strategy (pastastyle or kmh)",
                        required=False, default="pastastyle")
    
    parser.add_argument("--decompskeletonsize", type=int,
                        help="Number of skeleton sequences for the initial decomposition strategy",
                        required=False, default=300)
    
    parser.add_argument("--datatype", type=str,
                        help="Data type (dna, rna, or protein). Will be inferred if not provided",
                        required=False, default=None)
    
    parser.add_argument("--graphbuildmethod", type=str,
                        help="Method for building the alignment graph (mafft, mafftmerge, or initial)",
                        required=False, default="mafft")
    
    parser.add_argument("--graphbuildrestrict", type=str,
                        help="Prevent the alignment graph from adding edges that violate subalignments (true or false)",
                        required=False, default="False")
    
    parser.add_argument("--graphbuildhmmextend", type=str,
                        help="Extend the alignment graph MAFFT backbones with hmmer (true or false)",
                        required=False, default="False")
    
    parser.add_argument("--graphclustermethod", type=str,
                        help="Method for initial clustering of the alignment graph (mcl or none)",
                        required=False, default="mcl")
    
    parser.add_argument("--graphtracemethod", type=str,
                        help="Method for finding a trace from the alignment graph (minclusters, fm, mwtgreedy, or mwtsearch)",
                        required=False, default="minclusters")
    
    parser.add_argument("--graphtraceoptimize", type=str,
                        help="Run an optimization step on the graph trace (true or false)",
                        required=False, default="False")
    
    parser.add_argument("-r", "--mafftruns", type=int,
                        help="Number of MAFFT runs", required=False, default=10)
    
    parser.add_argument("-m", "--mafftsize", type=int,
                        help="Maximum size of MAFFT alignments", required=False, default=200)
    
    parser.add_argument("-f", "--inflationfactor", type=float,
                        help="MCL inflation factor", required=False, default=4)
    
    parser.add_argument("-c", "--constrain", type=str,
                        help="Constrain MAGUS to respect subalignments (true or false)", required=False, default="true")
    
    parser.add_argument("--onlyguidetree", type=str,
                        help="Only output the guide tree (true or false)", required=False, default="false")
    
    parser.add_argument("--recurse", type=str,
                        help="Allow MAGUS to recurse on large subsets (true or false)", required=False, default="true")
    
    parser.add_argument("--recurseguidetree", type=str,
                        help="If recursing, passes this argument as the guide tree option to the lower levels. (Default fasttree)", required=False, default="fasttree")
    
    parser.add_argument("--recursethreshold", type=int,
                        help="MAGUS will recursively align subsets above this threshold size", required=False, default=200)
    
    parser.add_argument("--alignsizelimit", type=float,
                        help="Size threshold for alignment compression (in GB)", required=False, default=100)
       
    return parser.parse_args()

if __name__ == '__main__':
    main()