'''
Created on Sep 29, 2020

@author: Vlad
'''

import os
import importlib
import json
import traceback
from magus_configuration import Configs
from magus_tasks import controller

'''
Tasks are self-contained, parallelizable units of work. 
For example, running MAFFT, compressing a subalignment, etc..
Primarily used to thread-parallelize MAFFT runs and node-parallelize subalignment operations.
Saved as JSON in task files, which are then read back by computing nodes with available threads.
This also serves the purpose of allowing aborted MAGUS runs to pick up where they left off.
'''

class Task:
    
    functionModuleMap = {"runCommand" : "magus_tools.external_tools",
                         "runAlignmentTask" : "magus_align.aligner",
                         "recordGapCounts" : "magus_align.merge.alignment_writer",
                         "buildInducedSubalignment" : "magus_align.merge.alignment_writer",
                         "compressSubalignment" : "magus_align.merge.alignment_writer"}
    
    def __init__(self, taskType, outputFile, taskArgs, **kwargs):
        self.taskType = taskType
        self.outputFile = outputFile
        self.taskArgs = taskArgs
        
        for attr in kwargs:
            vars(self)[attr] = kwargs.get(attr)    
        self.attributes =  list(vars(self).keys())
        
        self.isFinished = False
        self.future = None
        self.json = self.toJson()   
    
    def submitTask(self):
        submitTasks([self])
    
    def awaitTask(self):
        awaitTasks([self])
        
    def submitAndAwaitTask(self):
        self.submitTask()
        self.awaitTask()
        
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
            Configs.error("Task for {} threw an exception:\n{}".format(self.outputFile, exc))
            Configs.error(traceback.format_exc())
            raise
        finally:
            self.isFinished = True
        
    def checkFinished(self):
        if not self.isFinished:
            return False
        if self.future is not None:
            self.future.result()
        return True
        
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
    