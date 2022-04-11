'''
Created on Apr 14, 2020

@author: Vlad
'''

from gettext import find
import os
import time
from shutil import which
from sys import platform
import stat
from magus_helpers import sequenceutils
from os.path import basename
from glob import glob

def retrieve_packaged_binary(p):
    if platform == "linux" or platform == "linux2":
        for executable in glob(os.path.dirname(p) + "/**/*",recursive=True):
            if not os.path.isfile(executable):
                continue
            os.chmod(executable, os.stat(p).st_mode | stat.S_IEXEC)
        return p
    else:
        return None

def find_binary(rel_default_path):
    default_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), rel_default_path)
    return which(basename(rel_default_path)) or retrieve_packaged_binary(default_path)


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
    
    clustalPath = find_binary("magus_tools/clustal/clustalo")
    mafftPath = find_binary("magus_tools/mafft/mafft")
    mclPath = find_binary("magus_tools/mcl/bin/mcl")
    mlrmclPath = find_binary("magus_tools/mlrmcl/mlrmcl")
    hmmalignPath = find_binary("magus_tools/hmmer/hmmalign")
    hmmbuildPath = find_binary("magus_tools/hmmer/hmmbuild")
    hmmsearchPath = find_binary("magus_tools/hmmer/hmmsearch")
    fasttreePath = find_binary("magus_tools/fasttree/FastTreeMP")
    raxmlPath = find_binary("magus_tools/raxmlng/raxml-ng")
    
    logPath = None
    errorPath = None
    debugPath = None
    
    numCores = 1
    searchHeapLimit = 5000
    alignmentSizeLimit = 100
    
    @staticmethod
    def log(msg, path = None):
        print(msg)
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
    Configs.outputPath = os.path.abspath(args.output)
    
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
    
    Configs.logPath = os.path.join(Configs.workingDir, "log.txt")    
    Configs.errorPath = os.path.join(Configs.workingDir, "log_errors.txt")
    Configs.debugPath = os.path.join(Configs.workingDir, "log_debug.txt")
    
    Configs.alignmentSizeLimit = args.alignsizelimit
