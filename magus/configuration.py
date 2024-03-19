'''
Created on Apr 14, 2020

@author: Vlad
'''

import os
import time
from sys import stderr

from .helpers import sequenceutils


def find_binary(rel_path):
    import magus.tools
    default_path = os.path.join(os.path.dirname(os.path.abspath(magus.tools.__file__)), rel_path)
    return default_path


class Configs:
    
    workingDir = None
    sequencesPath = None
    subsetPaths = None
    subalignmentPaths = None
    backbonePaths = None
    guideTree = "fasttree"
    outputPath = None
    dataType = None
    
    decompositionMaxNumSubsets = 25
    decompositionMaxSubsetSize = 50
    decompositionStrategy = "pastastyle"
    decompositionSkeletonSize = 300
    #decompositionKmhIterations = 1
    
    graphBuildMethod = "mafft"
    graphBuildHmmExtend = False
    graphBuildRestrict = False
    graphClusterMethod = "mcl" 
    graphTraceMethod = "minclusters"
    graphTraceOptimize = False
    
    mafftRuns = 10
    mafftSize = 200
    mclInflationFactor = 4
    
    constrain = True
    onlyGuideTree = False
    recurse = True
    recurseGuideTree = "fasttree"
    recurseThreshold = 200
    
    clustalPath = find_binary("clustal/clustalo")
    mafftPath = find_binary("mafft/mafft")
    mclPath = find_binary("mcl/mcl")
    mlrmclPath = find_binary("mlrmcl/mlrmcl")
    hmmalignPath = find_binary("hmmer/hmmalign")
    hmmbuildPath = find_binary("hmmer/hmmbuild")
    hmmsearchPath = find_binary("hmmer/hmmsearch")
    fasttreePath = find_binary("fasttree/FastTreeMP")
    raxmlPath = find_binary("raxmlng/raxml-ng")

    overwrite = False
    
    logPath = None
    errorPath = None
    debugPath = None
    
    numCores = 1
    searchHeapLimit = 5000
    alignmentSizeLimit = 100
    
    @staticmethod
    def log(msg, path = None):
        print(msg, file=stderr)
        path = Configs.logPath if path is None else path
        Configs.writeMsg(msg, path)
    
    @staticmethod
    def error(msg, path = None):
        Configs.log(msg)
        path = Configs.errorPath if path is None else path
        Configs.writeMsg(msg, path)
    
    @staticmethod
    def debug(msg, path = None):
        path = Configs.debugPath if path is None else path
        Configs.writeMsg(msg, path)
    
    @staticmethod
    def writeMsg(msg, path):
        if path is not None:
            with open(path, 'a') as logFile:
                logFile.write("{}    {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S"), msg))
    
    @staticmethod
    def inferDataType(sequencesFile):
        if Configs.dataType is None:
            Configs.dataType = sequenceutils.inferDataType(sequencesFile)
            Configs.log("Data type wasn't specified. Inferred data type {} from {}".format(Configs.dataType.upper(), sequencesFile))
        return Configs.dataType 

def buildConfigs(args):
    Configs.outputPath = os.path.abspath(args.output) if args.output != '-' else "/dev/stdout"
    
    if args.directory is not None:
        Configs.workingDir = os.path.abspath(args.directory) 
    else:
        Configs.workingDir = os.path.join(os.path.dirname(Configs.outputPath), "magus_working_dir")
    if not os.path.exists(Configs.workingDir):
        os.makedirs(Configs.workingDir)
    
    Configs.sequencesPath = os.path.abspath(args.sequences) if args.sequences is not None else Configs.sequencesPath
    
    Configs.guideTree = os.path.abspath(args.guidetree) if args.guidetree is not None else Configs.guideTree
    if args.guidetree is not None:
        Configs.guideTree = os.path.abspath(args.guidetree) if os.path.exists(os.path.abspath(args.guidetree)) else args.guidetree
    
    Configs.subalignmentPaths = []
    for p in args.subalignments:
        path = os.path.abspath(p)
        if os.path.isdir(path):
            for filename in os.listdir(path):
                Configs.subalignmentPaths.append(os.path.join(path, filename))
        else:
            Configs.subalignmentPaths.append(path)
    
    Configs.backbonePaths = []
    for p in args.backbones:
        path = os.path.abspath(p)
        if os.path.isdir(path):
            for filename in os.listdir(path):
                Configs.backbonePaths.append(os.path.join(path, filename))
        else:
            Configs.backbonePaths.append(path)

    if args.numprocs > 0:
        Configs.numCores = args.numprocs
    else:
        Configs.numCores = os.cpu_count()

    Configs.decompositionMaxSubsetSize = args.maxsubsetsize
    Configs.decompositionMaxNumSubsets = args.maxnumsubsets
    Configs.decompositionStrategy = args.decompstrategy
    Configs.decompositionSkeletonSize = args.decompskeletonsize
    Configs.dataType = args.datatype
    
    Configs.graphBuildMethod = args.graphbuildmethod
    Configs.graphBuildHmmExtend = args.graphbuildhmmextend.lower() == "true"
    Configs.graphBuildRestrict = args.graphbuildrestrict.lower() == "true"
    Configs.graphClusterMethod = args.graphclustermethod
    Configs.graphTraceMethod = args.graphtracemethod
    Configs.graphTraceOptimize = args.graphtraceoptimize.lower() == "true"

    Configs.mafftRuns = args.mafftruns
    Configs.mafftSize = args.mafftsize
    Configs.mclInflationFactor = args.inflationfactor
    
    Configs.constrain = args.constrain.lower() == "true"
    Configs.onlyGuideTree = args.onlyguidetree.lower() == "true"
    Configs.recurse = args.recurse.lower() == "true"
    Configs.recurseGuideTree = args.recurseguidetree
    Configs.recurseThreshold = args.recursethreshold

    Configs.overwrite = args.overwrite
    
    Configs.logPath = os.path.join(Configs.workingDir, "log.txt")    
    Configs.errorPath = os.path.join(Configs.workingDir, "log_errors.txt")
    Configs.debugPath = os.path.join(Configs.workingDir, "log_debug.txt")
    
    Configs.alignmentSizeLimit = args.alignsizelimit
