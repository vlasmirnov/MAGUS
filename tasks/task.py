'''
Created on Sep 29, 2020

@author: Vlad
'''

import os
import importlib
import json
import traceback
from configuration import Configs
from tasks import controller


class Task:
    
    functionModuleMap = {"runCommand" : "tools.external_tools",
                         "runAlignmentTask" : "align.aligner",
                         "recordGapCounts" : "align.merge.alignment_writer",
                         "buildInducedSubalignment" : "align.merge.alignment_writer"}
    
    def __init__(self, taskType, outputFile, taskArgs, **kwargs):
        self.taskType = taskType
        self.outputFile = outputFile
        self.taskArgs = taskArgs
        
        for attr in kwargs:
            vars(self)[attr] = kwargs.get(attr)    
        self.attributes =  list(vars(self).keys())
           
        self.isFinished = False
        self.future = None
    
    def submitTask(self):
        submitTasks([self])
    
    def awaitTask(self):
        awaitTasks([self])
        
    def run(self):
        try:
            if not os.path.exists(self.outputFile):
                Configs.log("Running a task, output file: {}".format(self.outputFile))
                mod = importlib.import_module(Task.functionModuleMap[self.taskType])
                func = getattr(mod, self.taskType)
                func(**self.taskArgs)
                Configs.log("Completed a task, output file: {}".format(self.outputFile))
            else:
                Configs.log("File already exists: {}".format(self.outputFile))
        except Exception as exc:
            Configs.log("Task for {} threw an exception:\n{}".format(self.outputFile, exc))
            Configs.log(traceback.format_exc())
            raise
        
    def checkFinished(self):
        if not self.isFinished:
            return False
        try:
            if self.future is not None:
                self.future.result()
            return True
        except Exception as exc:
            Configs.log("Worker for task {} threw an exception:\n{}".format(self.outputFile, exc))
            Configs.log(traceback.format_exc())
            raise
    
    def toJson(self):
        mapper = {attr : getattr(self, attr) for attr in self.attributes}
        return json.dumps(mapper)
    
    def __eq__(self, other):
        if isinstance(other, Task):
            return self.outputFile == other.outputFile
        return NotImplemented

    def __hash__(self):
        return hash(self.outputFile)


def asCompleted(tasks):
    yield from controller.asCompleted(tasks)
            
def awaitTasks(tasks):
    controller.awaitTasks(tasks)
    
def submitTasks(tasks):
    controller.submitTasks(tasks)
    
