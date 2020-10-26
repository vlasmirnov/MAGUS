'''
Created on May 29, 2020

@author: Vlad
'''

import os
import time
import json

from align.decompose.decomposer import decomposeSequences 
from align.merge.merger import mergeSubalignments
from tools import external_tools
from configuration import Configs
from helpers import tasks, sequenceutils


def mainAlignmentTask():    
    mainTask = AlignTask(workingDir = Configs.workingDir, outputFile = Configs.outputPath,
                                    subsetPaths = Configs.subsetPaths, sequencesPath = Configs.sequencesPath,
                                    backbonePaths = Configs.backbonePaths, guideTreePath = Configs.guideTreePath)
    
    AlignTask.submittedTasksFile = os.path.join(Configs.workingDir, "alignment_tasks_submitted.txt")
    AlignTask.takenTasksFile = os.path.join(Configs.workingDir, "alignment_tasks_taken.txt")
    AlignTask.lockTasksFile = os.path.join(Configs.workingDir, "alignment_tasks_files.lock")
    
    submitMagusTasks([mainTask])
    runOrAwaitMagusTasks([mainTask])

def runAlignmentTask(task):
    if not os.path.exists(task.workingDir):
        os.makedirs(task.workingDir)
    
    if task.sequencesPath is not None:
        Configs.log("Aligning sequences {}".format(task.sequencesPath))
    
    if task.subsetPaths is None or len(task.subsetPaths) == 0:    
        decomposeSequences(task)
        alignSubsets(task)
    else:
        task.subalignmentPaths = task.subsetPaths
           
    mergeSubalignments(task)
    
def alignSubsets(task):
    Configs.log("Building {} subalignments..".format(len(task.subsetPaths)))
    subalignDir = os.path.join(task.workingDir, "subalignments")
    if not os.path.exists(subalignDir):
        os.makedirs(subalignDir)
    
    for file in task.subsetPaths:
        subset = sequenceutils.readFromFasta(file)
        if len(subset) <= Configs.decompositionMaxSubsetSize:
            Configs.log("Subset has {}/{} sequences, aligning with MAFFT..".format(len(subset), Configs.decompositionMaxSubsetSize))
            subalignmentPath = os.path.join(subalignDir, "subalignment_{}".format(os.path.basename(file)))
            task.subalignmentMafftTasks.append(external_tools.buildMafftAlignment(file, subalignmentPath))
        else:
            Configs.log("Subset has {}/{} sequences, recursively subaligning with MAGUS..".format(len(subset), Configs.decompositionMaxSubsetSize))
            subalignmentPath = os.path.join(subalignDir, "subalignment_{}".format(os.path.basename(file))) 
            subalignmentDir = os.path.join(subalignDir, os.path.splitext(os.path.basename(subalignmentPath))[0])
            task.subalignmentMagusTasks.append(AlignTask(outputFile = subalignmentPath, workingDir = subalignmentDir, sequencesPath = file))   
        task.subalignmentPaths.append(subalignmentPath)
            
    tasks.submitTasks(task.subalignmentMafftTasks)
    submitMagusTasks(task.subalignmentMagusTasks)
    #runOrAwaitMagusTasks(task.subalignmentMagusTasks, task.subalignmentMafftTasks)
    
    Configs.log("Prepared {} MAFFT and {} MAGUS tasks..".format(len(task.subalignmentMafftTasks), len(task.subalignmentMagusTasks)))
    
def runOrAwaitMagusTasks(alignTasks, additionalTasks = []):    
    while True:
        unfinishedTasks = [task for task in alignTasks if not os.path.exists(task.outputFile)]
        if len(unfinishedTasks) == 0:
            unfinishedAdditionals = [task for task in additionalTasks if not os.path.exists(task.outputFile)]
            if len(unfinishedAdditionals) == 0:
                return
        
        runningTask = None
        
        try:
            tasks.lockFiles(AlignTask.lockTasksFile)
            takenTasks = tasks.readTakenTasksFromFile(AlignTask.takenTasksFile)
            for task in unfinishedTasks:
                if task.outputFile not in takenTasks:
                    runningTask = task
                    break
            if runningTask is None:
                submittedTasks = tasks.readSubmittedTasksFromFile(takenTasks, 1, AlignTask.submittedTasksFile, AlignTask)
                runningTask = submittedTasks[0] if len(submittedTasks) > 0 else None
            if runningTask is not None:
                tasks.writeTakenTasksToFile([runningTask], AlignTask.takenTasksFile)
        finally:
            tasks.unlockFiles(AlignTask.lockTasksFile)
        
        if runningTask is not None:
            runAlignmentTask(runningTask)
        else:
            time.sleep(5)

def submitMagusTasks(alignTasks):
    try:
        tasks.lockFiles(AlignTask.lockTasksFile)
        takenTasks = tasks.readTakenTasksFromFile(AlignTask.takenTasksFile)
        submitTasks = []
        for task in alignTasks:
            if os.path.exists(task.outputFile):
                Configs.log("MAGUS alignment already exists: {}".format(task.outputFile))
            elif task.outputFile in takenTasks:
                Configs.log("MAGUS alignment already being built: {}".format(task.outputFile))
            else:
                submitTasks.append(task)
        tasks.writeSubmittedTasksToFile(submitTasks, AlignTask.submittedTasksFile)
    finally:
        tasks.unlockFiles(AlignTask.lockTasksFile)
    
    
class AlignTask:
    
    submittedTasksFile = None
    takenTasksFile = None
    lockTasksFile = None
    
    def __init__(self, **kwargs):
        self.outputFile = None
        self.workingDir = None
        self.sequencesPath = None
        self.subsetPaths = None
        self.backbonePaths = None
        self.guideTreePath = None
        
        for attr in vars(self):
            vars(self)[attr] = kwargs.get(attr)
        self.attributes =  list(vars(self).keys())  
    
        self.subalignmentPaths = []
        self.subalignmentMafftTasks = []
        self.subalignmentMagusTasks = []
        self.backboneTasks = []
        self.graph = None
    
    def toJson(self):
        mapper = {attr : getattr(self, attr) for attr in self.attributes}
        return json.dumps(mapper)

    def awaitSubalignments(self):
        runOrAwaitMagusTasks(self.subalignmentMagusTasks, self.subalignmentMafftTasks)


