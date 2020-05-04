'''
Created on Apr 14, 2020

@author: Vlad
'''

import os
import time

class Configs:
    
    workingDir = None
    sequencesPath = None
    subsetPaths = None
    backbonePaths = None
    guideTreePath = None
    outputPath = None
    
    mafftPath = None
    mafftRuns = 10
    mafftSize = 200
    mclPath = None
    mclInflationFactor = 4
    
    logPath = None
    graphPath = None
    clusterPath = None
    
    numCores = 1
    searchHeapLimit = 5000
    
    @staticmethod
    def log(msg):
        print(msg)
        if Configs.logPath is not None:
            with open(Configs.logPath, 'a') as logFile:
                logFile.write("{}    {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S"), msg))
    

def buildConfigs(args):
    Configs.outputPath = os.path.abspath(args.output)
    Configs.workingDir = os.path.abspath(args.directory) if args.directory is not None else os.path.dirname(Configs.outputPath)
    if not os.path.exists(Configs.workingDir):
        os.makedirs(Configs.workingDir)
    
    #Configs.sequencesPath = os.path.abspath(args.sequences)    
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
                
    #Configs.guideTreePath = os.path.abspath(args.guidetree)
    
    Configs.mafftPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tools/mafft/mafft")
    Configs.mafftRuns = args.mafftruns if args.mafftruns is not None else Configs.mafftRuns
    Configs.mafftSize = args.mafftsize if args.mafftsize is not None else Configs.mafftSize
    
    Configs.mclPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tools/mcl/bin/mcl")
    Configs.mclInflationFactor = args.inflationfactor if args.inflationfactor is not None else Configs.mclInflationFactor
    
    Configs.logPath = os.path.join(Configs.workingDir, "log.txt")
    Configs.graphPath = os.path.join(Configs.workingDir, "graph.txt")
    Configs.clusterPath = os.path.join(Configs.workingDir, "clusters.txt")
    
    Configs.numCores = os.cpu_count()
