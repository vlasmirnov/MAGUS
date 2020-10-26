'''
Created on Sep 29, 2020

@author: Vlad
'''

import subprocess
import os
import shutil
import random
import time
import concurrent.futures
import threading
import json
import traceback
import sys
from configuration import Configs

class Task:
    
    def __init__(self, **kwargs):
        self.outputFile = None
        self.command = None
        self.workingDir = None
        self.fileCopyMap = None
        
        for attr in vars(self):
            vars(self)[attr] = kwargs.get(attr)    
        self.attributes =  list(vars(self).keys())
           
        self.doneSignal = threading.Event()
        self.future = None
    
    
    def runTask(self):
        try:
            self.run()      
        finally:
            with TaskManager.managerLock:
                TaskManager.finishedTasks.add(self)
                TaskManager.runningTasks.remove(self)
                self.doneSignal.set()
                TaskManager.managerSignal.set()
        
    def run(self):
        if not os.path.exists(self.outputFile):
            Configs.log("Running a task, command: {}".format(self.command))
            subprocess.run(self.command, shell=True, cwd = self.workingDir)
            for srcPath, destPath in self.fileCopyMap.items():
                #shutil.copyfile(srcPath, destPath)
                shutil.move(srcPath, destPath)
            Configs.log("Completed a task, output file: {}".format(self.outputFile))
        
        
    def submitTask(self):
        submitTasks([self])
    
    def waitForTask(self):
        self.doneSignal.wait()     
        #Configs.log("Done signal from task with output file: {}".format(self.outputFile))
        try:
            if self.future is not None:
                self.future.result()
        except Exception as exc:
            Configs.log("Worker for task {} threw an exception:\n{}".format(self.command, exc))
            #exc_type, exc_value, exc_tb = sys.exc_info()
            #Configs.log(traceback.print_exception(exc_type, exc_value, exc_tb))
            Configs.log(traceback.format_exc())
            raise
        return self   
    
    def toJson(self):
        #attributes = ["outputFile", "command", "workingDir", "fileCopyMap"]
        mapper = {attr : getattr(self, attr) for attr in self.attributes}
        return json.dumps(mapper)
              

class TaskManager(threading.Thread):
    
    submittedTasksFile = None
    takenTasksFile = None
    lockTasksFile = None
    
    managerPool = None
    managerFuture = None
    managerSignal = threading.Event()
    managerLock = threading.Lock()
    managerStopSignal = False
    
    submittedTasks = set()
    waitingTasks = set()
    runningTasks = set()
    finishedTasks = set()
    taskPool = None
    threadsUsed = 0
    lastFilesCheckTime = 0
    error = None
    
    
def startTaskManager():
    Configs.log("Starting up the task manager..")
    TaskManager.submittedTasksFile = os.path.join(Configs.workingDir, "tasks_submitted.txt")
    TaskManager.takenTasksFile = os.path.join(Configs.workingDir, "tasks_taken.txt")
    TaskManager.lockTasksFile = os.path.join(Configs.workingDir, "tasks_files.lock")
    
    TaskManager.managerPool = concurrent.futures.ThreadPoolExecutor(max_workers = 1)
    TaskManager.taskPool = concurrent.futures.ThreadPoolExecutor(max_workers = Configs.numCores)
    TaskManager.managerFuture = TaskManager.managerPool.submit(runTaskManager)
    Configs.log("Task manager is up..")

def stopTaskManager():
    Configs.log("Winding down the task manager..")
    with TaskManager.managerLock:
        TaskManager.managerStopSignal = True
        TaskManager.managerSignal.set()
    TaskManager.managerFuture.result()
    TaskManager.managerPool.shutdown()
    Configs.log("Task manager stopped..")

def checkTaskManager():
    if not TaskManager.managerFuture.running():
        Configs.log("Task manager is dead for some reason..")
        TaskManager.managerFuture.result()
        raise Exception("Task manager is dead for some reason..")

def runTaskManager():
    while not TaskManager.managerStopSignal:
        with TaskManager.managerLock:
            
            if len(TaskManager.finishedTasks) > 0 or len(TaskManager.submittedTasks) > 0:
                Configs.log("Task manager status: {} submitted, {} waiting, {}/{} running, {} finished..".format(
                    len(TaskManager.submittedTasks), len(TaskManager.waitingTasks), len(TaskManager.runningTasks), Configs.numCores, len(TaskManager.finishedTasks)))
            
            dealWithFinishedTasks()
            dealWithSubmittedTasks()
            dealWithWaitingTasks()
            
            TaskManager.managerSignal.clear()
        TaskManager.managerSignal.wait(5)
    
    Configs.log("Waiting for {} tasks to finish..".format(len(TaskManager.runningTasks)))    
    TaskManager.taskPool.shutdown()

def dealWithFinishedTasks():
    for task in TaskManager.finishedTasks:
        TaskManager.threadsUsed = TaskManager.threadsUsed - 1
    TaskManager.finishedTasks = set()
    
def dealWithSubmittedTasks():
    threadsAvailable = Configs.numCores - TaskManager.threadsUsed
    if len(TaskManager.submittedTasks) > 0 or threadsAvailable > 0:
        try:
            lockFiles(TaskManager.lockTasksFile)
            takenTasks = readTakenTasksFromFile(TaskManager.takenTasksFile)
            
            runningList = []
            waitingList = []
            submittedList = []
            
            for task in TaskManager.submittedTasks:
                if os.path.exists(task.outputFile):
                    Configs.log("File already exists: {}".format(task.outputFile))
                    task.doneSignal.set()
                elif task.outputFile in takenTasks: 
                    Configs.log("Task already running: {}".format(task.outputFile))
                    waitingList.append(task)
                else:
                    submittedList.append(task)
                    if threadsAvailable > len(runningList):
                        runningList.append(task)
                    else:
                        waitingList.append(task)
            
            if threadsAvailable > len(runningList):
                additionalTasks = readSubmittedTasksFromFile(takenTasks, threadsAvailable - len(runningList), TaskManager.submittedTasksFile, Task)
                if len(additionalTasks) > 0:
                    Configs.log("Pulled {} tasks from the submitted tasks file..".format(len(additionalTasks)))
                runningList.extend(additionalTasks)
                
            writeSubmittedTasksToFile(submittedList, TaskManager.submittedTasksFile)    
            writeTakenTasksToFile(runningList, TaskManager.takenTasksFile)
        finally:
            unlockFiles(TaskManager.lockTasksFile)
        
        if len(runningList) + len(waitingList) > 0:
            Configs.log("Launching {} tasks and deferring {} tasks..".format(len(runningList), len(waitingList)))
        
        TaskManager.submittedTasks = set()    
        TaskManager.waitingTasks.update(waitingList)
        
        for task in runningList:
            TaskManager.threadsUsed = TaskManager.threadsUsed + 1
            TaskManager.runningTasks.add(task)
            task.future = TaskManager.taskPool.submit(task.runTask)
        
def dealWithWaitingTasks():
    timeSinceFileCheck = time.time() - TaskManager.lastFilesCheckTime
    if timeSinceFileCheck >= 5: 
        for task in list(TaskManager.waitingTasks):
            if os.path.exists(task.outputFile):
                Configs.log("Detected task completion: {}".format(task.outputFile))
                task.doneSignal.set()
                TaskManager.waitingTasks.remove(task)
        TaskManager.lastFilesCheckTime = time.time()
        

def asCompleted(tasks):
    with concurrent.futures.ThreadPoolExecutor() as waitPool:
        futures = {waitPool.submit(t.waitForTask) : t for t in tasks}
        for future in concurrent.futures.as_completed(futures):
            try:
                yield future.result()
            except Exception as exc:
                Configs.log("Worker for task {} threw an exception:\n{}".format(futures[future].command, exc))
                raise

#def addCallback(task, fn):
#    task.future.add_done_callback(fn)
    
def awaitTasks(tasks):
    with concurrent.futures.ThreadPoolExecutor() as waitPool:
        futures = {waitPool.submit(t.waitForTask) : t for t in tasks}
        #concurrent.futures.wait(futures)   
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as exc:
                Configs.log("Worker for task {} threw an exception:\n{}".format(futures[future].command, exc))
                raise
    
def submitTasks(tasks):
    with TaskManager.managerLock:
        checkTaskManager()
        TaskManager.submittedTasks.update(tasks)
        TaskManager.managerSignal.set()
    
def writeSubmittedTasksToFile(tasks, submittedTasksFile):
    with open(submittedTasksFile, 'a') as file:
        for task in tasks:
            file.write(task.toJson() + "\n")

def readSubmittedTasksFromFile(takenTasks, numTasks, submittedTasksFile, taskType):
    tasks = []
    if os.path.exists(submittedTasksFile) and numTasks > 0:
        with open(submittedTasksFile) as file:
            for line in file:
                mapper = json.loads(line.strip())
                task = taskType(**mapper)
                if task.outputFile not in takenTasks and not os.path.exists(task.outputFile):
                    tasks.append(task)
                    if len(tasks) >= numTasks:
                        break
    return tasks

def writeTakenTasksToFile(tasks, takenTasksFile):
    with open(takenTasksFile, 'a') as file:
        for task in tasks:
            file.write(task.outputFile + "\n")
            
def readTakenTasksFromFile(takenTasksFile):
    if not os.path.exists(takenTasksFile):
        return set()
    else:
        with open(takenTasksFile) as file:
            lines = [line.strip() for line in file]
        return set(lines)         

def lockFiles(file):
    while True:
        try:
            lock = open(file, 'x')
            lock.close()
            return
        except:
            time.sleep(random.random() + 1)
    
def unlockFiles(file):
    os.remove(file)