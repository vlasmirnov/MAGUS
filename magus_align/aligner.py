'''
Created on May 29, 2020

@author: Vlad
'''

import os
import shutil

from magus_align.alignment_context import AlignmentContext
from magus_align.decompose.decomposer import decomposeSequences
from magus_align.merge.merger import mergeSubalignments
from magus_tools import external_tools
from magus_configuration import Configs
from magus_helpers import sequenceutils
from magus_tasks import task

'''
Alignments are treated as "tasks", units of work that are written out to task files and 
processed as threads and/or compute nodes become available. 
MAGUS tasks will recursively generate MAGUS tasks over large subsets and MAFFT tasks over smaller subsets.
'''

def mainAlignmentTask():    
    args = {"workingDir" : Configs.workingDir, "outputFile" : Configs.outputPath,
            "subalignmentPaths" : Configs.subalignmentPaths, "sequencesPath" : Configs.sequencesPath,
            "backbonePaths" : Configs.backbonePaths, "guideTree" : Configs.guideTree}
    task = createAlignmentTask(args)
    task.submitTask()
    task.awaitTask()
    
def createAlignmentTask(args):
    return task.Task(taskType = "runAlignmentTask", outputFile = args["outputFile"], taskArgs = args)

def runAlignmentTask(**kwargs):
    '''
    The standard MAGUS task: 
    decompose the data into subsets, align each subset, and merge the subalignments.
    '''
    
    with AlignmentContext(**kwargs) as context:
        if context.sequencesPath is not None:
            Configs.log("Aligning sequences {}".format(context.sequencesPath))
        
        decomposeSequences(context)
        if Configs.onlyGuideTree:
            Configs.log("Outputting only the guide tree, as requested..")
            shutil.copyfile(os.path.join(context.workingDir, "decomposition", "initial_tree", "initial_tree.tre"), context.outputFile)
            return
        
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
        
    mafftThreshold = max(Configs.mafftSize, Configs.decompositionMaxSubsetSize, Configs.recurseThreshold)
    
    for file in context.subsetPaths:
        subset = sequenceutils.readFromFasta(file)
        subalignmentPath = os.path.join(subalignDir, "subalignment_{}".format(os.path.basename(file)))
        context.subalignmentPaths.append(subalignmentPath)
        
        if os.path.exists(subalignmentPath):
            Configs.log("Existing subalignment file detected: {}".format(subalignmentPath))       
             
        elif len(subset) <= mafftThreshold or not Configs.recurse:
            Configs.log("Subset has {}/{} sequences, aligning with MAFFT..".format(len(subset), mafftThreshold))            
            subalignmentTask = external_tools.buildMafftAlignment(file, subalignmentPath)
            context.subalignmentTasks.append(subalignmentTask)
            
        else:
            Configs.log("Subset has {}/{} sequences, recursively subaligning with MAGUS..".format(len(subset), mafftThreshold))
            subalignmentDir = os.path.join(subalignDir, os.path.splitext(os.path.basename(subalignmentPath))[0])
            subalignmentTask = createAlignmentTask({"outputFile" : subalignmentPath, "workingDir" : subalignmentDir, 
                                                    "sequencesPath" : file, "guideTree" : Configs.recurseGuideTree})   
            context.subalignmentTasks.append(subalignmentTask)
                
    task.submitTasks(context.subalignmentTasks)
    Configs.log("Prepared {} subset alignment tasks..".format(len(context.subalignmentTasks)))

