'''
Created on Apr 14, 2020

@author: Vlad
'''

import time
import argparse
import sys

from align.aligner import alignSequences
from align.merge.merger import mergeSubalignments
from configuration import buildConfigs, Configs


def main(args):   
    startTime = time.time()
    buildConfigs(args)
    Configs.log("GCM was run with: {}".format(" ".join(sys.argv)))
    
    if Configs.sequencesPath is not None:
        Configs.log("Aligning sequences {}".format(Configs.sequencesPath))
        alignSequences(Configs.workingDir, Configs.sequencesPath, Configs.outputPath)
    else:
        Configs.log("Merging {} sequences..".format(len(Configs.subsetPaths)))
        mergeSubalignments(Configs.workingDir, Configs.subsetPaths, Configs.outputPath)
    
    endTime = time.time()
    Configs.log("GCM finished in {} seconds..".format(endTime-startTime))
    

if __name__ == '__main__':
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
                        help="User guide tree for alignment",
                        required=False, default=None)
    
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
    
    parser.add_argument("-r", "--mafftruns", type=int,
                        help="Number of MAFFT runs", required=False, default=10)
    
    parser.add_argument("-m", "--mafftsize", type=int,
                        help="Maximum size of MAFFT alignments", required=False, default=200)
    
    parser.add_argument("-f", "--inflationfactor", type=float,
                        help="MCL inflation factor", required=False, default=4)

    main(parser.parse_args())