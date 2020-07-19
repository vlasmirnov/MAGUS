'''
Created on Apr 14, 2020

@author: Vlad
'''

import os
import time
from helpers import sequenceutils

class Configs:
    
    workingDir = None
    sequencesPath = None
    subsetPaths = None
    backbonePaths = None
    guideTreePath = None
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
    
    mafftRuns = 10
    mafftSize = 200
    mclInflationFactor = 4
    
    mafftPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tools/mafft/mafft")
    mclPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tools/mcl/bin/mcl")
    hmmalignPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tools/hmmer/hmmalign")
    hmmbuildPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tools/hmmer/hmmbuild")
    hmmsearchPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tools/hmmer/hmmsearch")
    fasttreePath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tools/fasttree/FastTree")
    raxmlPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tools/raxmlng/raxml-ng")
    
    logPath = None
    
    numCores = 1
    searchHeapLimit = 5000
    
    @staticmethod
    def log(msg):
        print(msg)
        if Configs.logPath is not None:
            with open(Configs.logPath, 'a') as logFile:
                logFile.write("{}    {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S"), msg))
    
    @staticmethod
    def inferDataType(sequencesFile):
        if Configs.dataType is None:
            Configs.dataType = sequenceutils.inferDataType(sequencesFile)
            Configs.log("Data type wasn't specified. Inferred data type {} from {}".format(Configs.dataType.upper(), sequencesFile))
        return Configs.dataType 

def buildConfigs(args):
    Configs.outputPath = os.path.abspath(args.output)
    
    if args.directory is not None:
        Configs.workingDir = os.path.abspath(args.directory) 
    else:
        Configs.workingDir = os.path.join(os.path.dirname(Configs.outputPath), "gcm_working_dir")
    if not os.path.exists(Configs.workingDir):
        os.makedirs(Configs.workingDir)
    
    Configs.sequencesPath = os.path.abspath(args.sequences) if args.sequences is not None else Configs.sequencesPath
    Configs.guideTreePath = os.path.abspath(args.guidetree) if args.guidetree is not None else Configs.guideTreePath
    
    Configs.subsetPaths = []
    for p in args.subalignments:
        path = os.path.abspath(p)
        if os.path.isdir(path):
            for filename in os.listdir(path):
                Configs.subsetPaths.append(os.path.join(path, filename))
        else:
            Configs.subsetPaths.append(path)
    
    Configs.backbonePaths = []
    for p in args.backbones:
        path = os.path.abspath(p)
        if os.path.isdir(path):
            for filename in os.listdir(path):
                Configs.backbonePaths.append(os.path.join(path, filename))
        else:
            Configs.backbonePaths.append(path)
    
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

    Configs.mafftRuns = args.mafftruns
    Configs.mafftSize = args.mafftsize
    Configs.mclInflationFactor = args.inflationfactor
    
    Configs.logPath = os.path.join(Configs.workingDir, "log.txt")    
    Configs.numCores = os.cpu_count()

       
