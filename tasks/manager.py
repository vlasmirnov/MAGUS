'''
Created on Oct 26, 2020

@author: Vlad
'''

import threading
import os
import time
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
    
    waitingTasks = {}
    submittedTasks = set()
    runningTasks = set()
    finishedTasks = set()
    failedTasks = set()
    
    taskPool = None
    threadsUsed = 0
    lastFilesCheckTime = 0
    lastDebugTime = 0
    serialTaskTypes = {"runAlignmentTask", "buildInducedSubalignment"}
    
    
def startTaskManager():
    Configs.log("Starting up the task manager..")
    
    tasksDir = os.path.join(Configs.workingDir, "tasks")
    if not os.path.exists(tasksDir):
        os.makedirs(tasksDir)  
    TaskManager.pendingTasksFile = os.path.join(tasksDir, "tasks_pending.txt")
    TaskManager.runningTasksFile = os.path.join(tasksDir, "tasks_running.txt")
    TaskManager.lockTasksFile = os.path.join(tasksDir, "tasks.lock")
    
    TaskManager.managerPool = concurrent.futures.ThreadPoolExecutor(max_workers = 1)
    TaskManager.managerFuture = TaskManager.managerPool.submit(runTaskManager)
    TaskManager.taskPool = concurrent.futures.ThreadPoolExecutor(max_workers = Configs.numCores)
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
        dealWithFinishedTasks()        
        TaskManager.managerPool.shutdown()
        Configs.log("Task manager stopped..")

def runTaskManager():
    try:
        while not TaskManager.managerStopSignal:
            with TaskManager.managerLock:
                dealWithErrors()
                dealWithFinishedTasks()
                dealWithPendingTasks()
                dealWithWaitingTasks()
                TaskManager.managerSignal.clear()
            TaskManager.managerSignal.wait(5)
    finally:
        TaskManager.observerSignal.set()  

def dealWithErrors():
    for task in TaskManager.failedTasks:
        Configs.error("Task manager found a failed task: {}".format(task.outputFile))
        if task.future is not None:
            task.future.result()
        
def dealWithFinishedTasks():
    stoppedRunning = set(t for t in TaskManager.finishedTasks | TaskManager.failedTasks if t.taskType == "runAlignmentTask")
    if len(stoppedRunning) > 0:
        with files.FileLock(TaskManager.lockTasksFile):
            runningTasks = files.readTasksFromFile(TaskManager.runningTasksFile)
            stillRunningTasks = [t for t in runningTasks if t not in stoppedRunning]
            if len(stillRunningTasks) < len(runningTasks):
                if len(stillRunningTasks) == 0 and os.path.exists(TaskManager.runningTasksFile):
                    os.remove(TaskManager.runningTasksFile)
                else:
                    files.writeTasksToFile(stillRunningTasks, TaskManager.runningTasksFile, append = False)
        
    if len(TaskManager.failedTasks) > 0:
        with files.FileLock(TaskManager.lockTasksFile):
            files.writeTasksToFile(TaskManager.failedTasks, TaskManager.pendingTasksFile, append = True)    
        
    TaskManager.finishedTasks = set()
    TaskManager.failedTasks = set()

def dealWithPendingTasks():
    numToLaunch = min(1, Configs.numCores - TaskManager.threadsUsed)
    newTasks = []
    for t in TaskManager.submittedTasks:
        if os.path.exists(t.outputFile):
            Configs.log("File already exists: {}".format(t.outputFile))
            t.isFinished = True
            TaskManager.observerSignal.set()
        else:
            newTasks.append(t)
            TaskManager.waitingTasks[t.outputFile] = t
    
    if len(newTasks) == 0 and (numToLaunch == 0 or not os.path.exists(TaskManager.pendingTasksFile)):
        return
    else:
        with files.FileLock(TaskManager.lockTasksFile):
            pendingTasks = files.readTasksFromFile(TaskManager.pendingTasksFile)
            runningTasks = files.readTasksFromFile(TaskManager.runningTasksFile)
            taskSet = set(pendingTasks) | set(runningTasks)
            newTasksToLaunch = newTasks
            if len(taskSet) > 0:
                newTasksToLaunch = [t for t in newTasks if t not in taskSet]
                
            launchedNewTasks, remainingNewTasks = launchTasks(newTasksToLaunch, numToLaunch)
            numToLaunch = numToLaunch - len(launchedNewTasks)
            launchedPendingTasks, remainingPendingTasks = launchTasks(pendingTasks, numToLaunch)
            
            writeRunningTasks = [t for t in launchedNewTasks + launchedPendingTasks if t.taskType == "runAlignmentTask"]
            if len(writeRunningTasks) > 0:
                files.writeTasksToFile(writeRunningTasks, TaskManager.runningTasksFile, append = True)
            
            if len(remainingPendingTasks) + len(remainingNewTasks) == 0 and os.path.exists(TaskManager.pendingTasksFile):
                os.remove(TaskManager.pendingTasksFile)
            else:
                if len(remainingPendingTasks) < len(pendingTasks):
                    files.writeTasksToFile(remainingPendingTasks, TaskManager.pendingTasksFile, append = False)
                files.writeTasksToFile(remainingNewTasks, TaskManager.pendingTasksFile, append = True)

        if len(launchedNewTasks) + len(launchedPendingTasks) + len(remainingNewTasks) + len(remainingPendingTasks) > 0:     
            Configs.log("Launched {} submitted tasks and {} pending tasks, deferred {} submitted tasks and {} pending tasks.."
                        .format(len(launchedNewTasks), len(launchedPendingTasks), len(remainingNewTasks), len(remainingPendingTasks)))
    
    TaskManager.submittedTasks = set()

def dealWithWaitingTasks():
    timeSinceFileCheck = time.time() - TaskManager.lastFilesCheckTime
    if timeSinceFileCheck >= 5: 
        for file, task in list(TaskManager.waitingTasks.items()):
            if os.path.exists(file):
                Configs.log("Detected task completion: {}".format(file))
                TaskManager.waitingTasks.pop(file)
                task.isFinished = True
                TaskManager.observerSignal.set()
        TaskManager.lastFilesCheckTime = time.time()
    
    timeSinceDebug = time.time() - TaskManager.lastDebugTime
    if timeSinceDebug >= 60:  
        TaskManager.lastDebugTime = time.time()
        for task in TaskManager.runningTasks:
            Configs.debug("Still running task {}, status {}".format(task.outputFile, task.future._state if task.future is not None else "N/A"))
        for file in TaskManager.waitingTasks:
            Configs.debug("Still waiting on task {}".format(file))

def launchTasks(tasks, numTasksToLaunch):
    launchedTasks = []
    remainingTasks = []    
    threadsAvailable = Configs.numCores - TaskManager.threadsUsed
    toLaunch = min(numTasksToLaunch, threadsAvailable)    
    for task in tasks:
        launched = False
        if toLaunch > 0 and task.taskType not in TaskManager.serialTaskTypes:
            task.future = TaskManager.taskPool.submit(runTask, task)
            launched = True
        elif toLaunch > 0 and TaskManager.observerWaiting and TaskManager.observerTask is None:
            TaskManager.observerTask = task
            TaskManager.observerSignal.set()
            launched = True
        
        if launched:
            Configs.log("Launched a new task.. {}/{} threads used, type: {}, output file: {}".format(TaskManager.threadsUsed, Configs.numCores, task.taskType, task.outputFile))
            toLaunch = toLaunch - 1
            launchedTasks.append(task)
        else:
            remainingTasks.append(task)
    return launchedTasks, remainingTasks

def runTask(task):
    with TaskManager.managerLock:
        if task.taskType not in TaskManager.serialTaskTypes:
            TaskManager.threadsUsed = TaskManager.threadsUsed + 1
        TaskManager.runningTasks.add(task)
    
    failed = False
    try:
        task.run()     
    except:
        failed = True
        raise
    finally:
        with TaskManager.managerLock:
            if task.taskType not in TaskManager.serialTaskTypes:
                TaskManager.threadsUsed = TaskManager.threadsUsed - 1
            TaskManager.runningTasks.remove(task)
            TaskManager.failedTasks.add(task) if failed else TaskManager.finishedTasks.add(task)
            if task.outputFile in TaskManager.waitingTasks:
                t = TaskManager.waitingTasks.pop(task.outputFile)
                t.future = task.future
                t.isFinished = True
            TaskManager.managerSignal.set()
            TaskManager.observerSignal.set()
                
        