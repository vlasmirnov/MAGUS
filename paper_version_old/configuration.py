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
    guideTreePath = None
    outputPath = None
    
    mafftPath = None
    mafftRuns = 10
    mafftSize = 200
    mclPath = None
    
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
    #Configs.workingDir = args.directory if args.directory else os.getcwd() 
    Configs.workingDir = args.directory if args.directory is not None else os.path.dirname(args.output)
    Configs.sequencesPath = args.sequences
    Configs.subsetPaths = args.subalignments
    #Configs.guideTreePath = args.guidetree
    Configs.outputPath = args.output
    
    Configs.mafftPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools/mafft/mafft")
    Configs.mafftRuns = args.mafftruns if args.mafftruns is not None else Configs.mafftRuns
    Configs.mafftSize = args.mafftsize if args.mafftsize is not None else Configs.mafftSize
    Configs.mclPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools/mcl/bin/mcl")
    
    Configs.logPath = os.path.join(Configs.workingDir, "log.txt")
    Configs.graphPath = os.path.join(Configs.workingDir, "graph.txt")
    Configs.clusterPath = os.path.join(Configs.workingDir, "clusters.txt")
    
    Configs.numCores = os.cpu_count()
