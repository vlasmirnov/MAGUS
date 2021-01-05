'''
Created on Oct 26, 2020

@author: Vlad
'''

import threading
import os
import time
import traceback
import concurrent.futures
from configuration import Configs
from tasks import files

class TaskManager():
    
    pendingTasksFile = None
    runningTasksFile = None
    lockTasksFile = None
    
    managerPool = None
    managerFuture = None
    managerSignal = threading.Event()
    managerLock = threading.Lock()
    managerStopSignal = False
    
    observerSignal = threading.Event()
    observerWaiting = False
    observerTask = None 
    
    submittedTasks = set()
    waitingTasks = set()
    runningTasks = set()
    finishedTasks = set()
    taskPool = None
    threadsUsed = 0
    futuresTaskMap = {}
    lastFilesCheckTime = 0
    serialTaskTypes = {"runAlignmentTask", "buildInducedSubalignment"}
    
    
def startTaskManager():
    Configs.log("Starting up the task manager..")
    TaskManager.pendingTasksFile = os.path.join(Configs.workingDir, "tasks_pending.txt")
    TaskManager.runningTasksFile = os.path.join(Configs.workingDir, "tasks_running.txt")
    TaskManager.lockTasksFile = os.path.join(Configs.workingDir, "tasks.lock")
    
    TaskManager.managerPool = concurrent.futures.ThreadPoolExecutor(max_workers = 1)
    TaskManager.managerFuture = TaskManager.managerPool.submit(runTaskManager)
    TaskManager.taskPool = concurrent.futures.ProcessPoolExecutor(max_workers = Configs.numCores)
    Configs.log("Task manager is up..")

def stopTaskManager():
    TaskManager.managerStopSignal = True
    with TaskManager.managerLock:
        TaskManager.managerSignal.set()
    try:
        Configs.log("Winding down the task manager..")
        TaskManager.managerFuture.result()
    finally:
        Configs.log("Waiting for {} tasks to finish..".format(len(TaskManager.runningTasks)))
        TaskManager.taskPool.shutdown()
        
        with files.FileLock(TaskManager.lockTasksFile):
            runningTasks = files.readTasksFromFile(TaskManager.runningTasksFile)
            checkRunningTasks(runningTasks)
        
        TaskManager.managerPool.shutdown()
        Configs.log("Task manager stopped..")

def runTaskManager():
    try:
        while True:
            with TaskManager.managerLock:
                if TaskManager.managerStopSignal:
                    break
                        
                if len(TaskManager.finishedTasks) > 0 or len(TaskManager.submittedTasks) > 0:
                    Configs.log("Task manager status: {} submitted, {} waiting, {}/{} running, {} finished..".format(
                        len(TaskManager.submittedTasks), len(TaskManager.waitingTasks), len(TaskManager.runningTasks), Configs.numCores, len(TaskManager.finishedTasks)))
                                
                dealWithPendingTasks()
                checkWaitingTasks()
                
                TaskManager.managerSignal.clear()
            TaskManager.managerSignal.wait(5)
    finally:
        TaskManager.observerSignal.set()  
        
def taskCallback(future):
    task = TaskManager.futuresTaskMap.pop(future)
    try:
        task.future = future
        with TaskManager.managerLock:
            TaskManager.finishedTasks.add(task)
            TaskManager.runningTasks.remove(task)
            TaskManager.threadsUsed = TaskManager.threadsUsed - 1
            TaskManager.managerSignal.set()
        setTaskFinished(task)
    except Exception as exc:
        Configs.log("Task callback threw an exception:\n{}".format(exc))
        Configs.log(traceback.format_exc())
        Configs.log("Task file: {}".format(task.outputFile))
        raise

def setTaskFinished(task):
    task.isFinished = True
    TaskManager.observerSignal.set()

def dealWithPendingTasks():
    threadsAvailable = Configs.numCores - TaskManager.threadsUsed
    if len(TaskManager.submittedTasks) > 0 or threadsAvailable > 0:
        
        with files.FileLock(TaskManager.lockTasksFile):
            runningTasks = files.readTasksFromFile(TaskManager.runningTasksFile)
            pendingTasks = files.readTasksFromFile(TaskManager.pendingTasksFile)
            newTasks = []
            
            for task in TaskManager.submittedTasks:
                if os.path.exists(task.outputFile):
                    Configs.log("File already exists: {}".format(task.outputFile))
                    setTaskFinished(task)
                elif task in pendingTasks or task in runningTasks: 
                    Configs.log("Task already pending or running: {}".format(task.outputFile))
                    TaskManager.waitingTasks.add(task)
                else:
                    newTasks.append(task)
            
            checkRunningTasks(runningTasks)
            newRunningTasks, newPendingTasks = assignTasks(newTasks, pendingTasks)
        
        if len(newRunningTasks) + len(newPendingTasks) > 0:
            Configs.log("Launching {} tasks and deferring {} tasks..".format(len(newRunningTasks), len(newPendingTasks)))
        
        TaskManager.finishedTasks = set()
        TaskManager.submittedTasks = set()    
        
        for task in newRunningTasks: 
            Configs.log("New running task, type {}, output file {}".format(task.taskType, task.outputFile))
            if task.taskType not in TaskManager.serialTaskTypes:
                TaskManager.threadsUsed = TaskManager.threadsUsed + 1
                TaskManager.runningTasks.add(task)
                future = TaskManager.taskPool.submit(task.run)
                Configs.log("Added task to task pool ({} threads used), type {}, output file {}".format(TaskManager.threadsUsed, task.taskType, task.outputFile))
                TaskManager.futuresTaskMap[future] = task
                future.add_done_callback(taskCallback)
                Configs.log("Future status {} for file {}".format(future._state, task.outputFile))
        for task in newPendingTasks:
            TaskManager.waitingTasks.add(task)

def assignTasks(newTasks, pendingTasks):
    stillPendingTasks, newRunningTasks, newPendingTasks = [], [], []
    
    threadsAvailable = Configs.numCores - TaskManager.threadsUsed
    for num, task in enumerate(newTasks + pendingTasks):
        if threadsAvailable > 0 and task.taskType not in TaskManager.serialTaskTypes:
            newRunningTasks.append(task)
            threadsAvailable = threadsAvailable - 1
        elif task.taskType in TaskManager.serialTaskTypes and TaskManager.observerWaiting and TaskManager.observerTask is None:
            newRunningTasks.append(task)
            TaskManager.observerTask = task
            TaskManager.observerSignal.set()
        elif num < len(newTasks):
            newPendingTasks.append(task) 
        else:
            stillPendingTasks.append(task)
    
    if len(stillPendingTasks) < len(pendingTasks):
        files.writeTasksToFile(stillPendingTasks, TaskManager.pendingTasksFile, append = False)
    files.writeTasksToFile(newRunningTasks, TaskManager.runningTasksFile, append = True)
    files.writeTasksToFile(newPendingTasks, TaskManager.pendingTasksFile, append = True)
    return newRunningTasks, newPendingTasks

def checkRunningTasks(runningTasks):
    stillRunningTasks = [t for t in runningTasks if t not in TaskManager.finishedTasks]
    if len(stillRunningTasks) < len(runningTasks):
        files.writeTasksToFile(stillRunningTasks, TaskManager.runningTasksFile, append = False)
        
def checkWaitingTasks():
    timeSinceFileCheck = time.time() - TaskManager.lastFilesCheckTime
    if timeSinceFileCheck >= 5: 
        for task in list(TaskManager.waitingTasks):
            if os.path.exists(task.outputFile):
                Configs.log("Detected task completion: {}".format(task.outputFile))
                TaskManager.waitingTasks.remove(task)
                setTaskFinished(task)
        TaskManager.lastFilesCheckTime = time.time()
        
