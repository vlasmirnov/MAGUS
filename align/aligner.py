'''
Created on May 29, 2020

@author: Vlad
'''

import os

from align.alignment_context import AlignmentContext
from align.decompose.decomposer import decomposeSequences 
from align.merge.merger import mergeSubalignments
from tools import external_tools
from configuration import Configs
from helpers import sequenceutils
from tasks import task


def mainAlignmentTask():    
    args = {"workingDir" : Configs.workingDir, "outputFile" : Configs.outputPath,
            "subalignmentPaths" : Configs.subalignmentPaths, "sequencesPath" : Configs.sequencesPath,
            "backbonePaths" : Configs.backbonePaths, "guideTreePath" : Configs.guideTreePath}
    task = createAlignmentTask(args)
    task.submitTask()
    task.awaitTask()
    
def createAlignmentTask(args):
    return task.Task(taskType = "runAlignmentTask", outputFile = args["outputFile"], taskArgs = args)

def runAlignmentTask(**kwargs):
    context = AlignmentContext(**kwargs)

    if context.sequencesPath is not None:
        Configs.log("Aligning sequences {}".format(context.sequencesPath))
    
    decomposeSequences(context)
    alignSubsets(context)
    mergeSubalignments(context)
    
def alignSubsets(context):
    if len(context.subalignmentPaths) > 0:
        Configs.log("Subalignment paths already provided, skipping subalignments..")
        return
    
    Configs.log("Building {} subalignments..".format(len(context.subsetPaths)))
    subalignDir = os.path.join(context.workingDir, "subalignments")
    if not os.path.exists(subalignDir):
        os.makedirs(subalignDir)
    
    for file in context.subsetPaths:
        subset = sequenceutils.readFromFasta(file)
        subalignmentPath = os.path.join(subalignDir, "subalignment_{}".format(os.path.basename(file)))
        context.subalignmentPaths.append(subalignmentPath)
        
        if os.path.exists(subalignmentPath):
            Configs.log("Existing subalignment file detected: {}".format(subalignmentPath))       
             
        elif len(subset) <= Configs.mafftSize: #Configs.decompositionMaxSubsetSize:
            Configs.log("Subset has {}/{} sequences, aligning with MAFFT..".format(len(subset), Configs.mafftSize))            
            subalignmentTask = external_tools.buildMafftAlignment(file, subalignmentPath)
            context.subalignmentTasks.append(subalignmentTask)
            
        else:
            Configs.log("Subset has {}/{} sequences, recursively subaligning with MAGUS..".format(len(subset), Configs.mafftSize))
            subalignmentDir = os.path.join(subalignDir, os.path.splitext(os.path.basename(subalignmentPath))[0])
            subalignmentTask = createAlignmentTask({"outputFile" : subalignmentPath, "workingDir" : subalignmentDir, "sequencesPath" : file})   
            context.subalignmentTasks.append(subalignmentTask)
                
    task.submitTasks(context.subalignmentTasks)
    Configs.log("Prepared {} subset alignment tasks..".format(len(context.subalignmentTasks)))

