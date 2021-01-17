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
    
    submittedTasks = set()
    waitingTasks = set()
    runningTasks = set()
    finishedTasks = set()
    failedTasks = set()
    
    taskPool = None
    threadsUsed = 0
    futuresTaskMap = {}
    lastFilesCheckTime = 0
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
        dealWithFinishedTasks()        
        TaskManager.managerPool.shutdown()
        Configs.log("Task manager stopped..")

def runTaskManager():
    try:
        while not TaskManager.managerStopSignal:
            with TaskManager.managerLock:
                dealWithFinishedTasks()
                dealWithPendingTasks()
                dealWithWaitingTasks()
                TaskManager.managerSignal.clear()
            TaskManager.managerSignal.wait(5)
    finally:
        TaskManager.observerSignal.set()  
        
def dealWithFinishedTasks():
    doneTasks = [t for t in TaskManager.finishedTasks if t.taskType == "runAlignmentTask"]
    failedTasks = [t for t in TaskManager.failedTasks if t.taskType == "runAlignmentTask"]
    TaskManager.finishedTasks = set()
    TaskManager.failedTasks = set()
    
    if len(doneTasks) + len(failedTasks) > 0:
        with files.FileLock(TaskManager.lockTasksFile):
            runningTasks = files.readTasksFromFile(TaskManager.runningTasksFile)
            taskSet = set(doneTasks) | set(failedTasks)
            stillRunningTasks = [t for t in runningTasks if t not in taskSet]
            
            if len(stillRunningTasks) < len(runningTasks):
                if len(stillRunningTasks) == 0 and os.path.exists(TaskManager.runningTasksFile):
                    os.remove(TaskManager.runningTasksFile)
                else:
                    files.writeTasksToFile(stillRunningTasks, TaskManager.runningTasksFile, append = False)
                    
            if len(failedTasks) > 0:
                files.writeTasksToFile(failedTasks, TaskManager.pendingTasksFile, append = True)    

def dealWithPendingTasks():
    numToLaunch = min(1, Configs.numCores - TaskManager.threadsUsed)
    newTasks = []
    for t in TaskManager.submittedTasks:
        if os.path.exists(t.outputFile):
            Configs.log("File already exists: {}".format(t.outputFile))
            setTaskFinished(t)
        else:
            newTasks.append(t)
    
    if len(newTasks) == 0 and (numToLaunch == 0 or not os.path.exists(TaskManager.pendingTasksFile)):
        return
    else:
        with files.FileLock(TaskManager.lockTasksFile):
            pendingTasks = files.readTasksFromFile(TaskManager.pendingTasksFile)
            runningTasks = files.readTasksFromFile(TaskManager.runningTasksFile)
            taskSet = set(pendingTasks) | set(runningTasks)
            newTasksToLaunch = newTasks
            if len(taskSet) > 0:
                newTasksToLaunch = []
                for t in newTasks:
                    if t in taskSet:
                        Configs.log("Task already pending or running: {}".format(t.outputFile))
                        TaskManager.waitingTasks.add(t)
                    else:
                        newTasksToLaunch.append(t)
                
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

        TaskManager.waitingTasks.update(remainingNewTasks)    
        if len(launchedNewTasks) + len(launchedPendingTasks) + len(remainingNewTasks) + len(remainingPendingTasks) > 0:     
            Configs.log("Launched {} submitted tasks and {} pending tasks, deferred {} submitted tasks and {} pending tasks.."
                        .format(len(launchedNewTasks), len(launchedPendingTasks), len(remainingNewTasks), len(remainingPendingTasks)))
    
    TaskManager.submittedTasks = set()

def dealWithWaitingTasks():
    timeSinceFileCheck = time.time() - TaskManager.lastFilesCheckTime
    if timeSinceFileCheck >= 5: 
        for task in list(TaskManager.waitingTasks):
            if os.path.exists(task.outputFile):
                Configs.log("Detected task completion: {}".format(task.outputFile))
                TaskManager.waitingTasks.remove(task)
                setTaskFinished(task)
        TaskManager.lastFilesCheckTime = time.time()

def launchTasks(tasks, numTasksToLaunch):
    launchedTasks = []
    remainingTasks = []    
    threadsAvailable = Configs.numCores - TaskManager.threadsUsed
    toLaunch = min(numTasksToLaunch, threadsAvailable)    
    for task in tasks:
        launchable = task.taskType not in TaskManager.serialTaskTypes or (TaskManager.observerWaiting and TaskManager.observerTask is None)
        if toLaunch > 0 and launchable:
            launchTask(task)
            toLaunch = toLaunch - 1
            launchedTasks.append(task)
        else:
            remainingTasks.append(task)
    return launchedTasks, remainingTasks

def launchTask(task):
    Configs.log("Launching a new task.. {}/{} threads used, type: {}, output file: {}".format(TaskManager.threadsUsed, Configs.numCores, task.taskType, task.outputFile))
    TaskManager.runningTasks.add(task)
    if task.taskType not in TaskManager.serialTaskTypes:
        TaskManager.threadsUsed = TaskManager.threadsUsed + 1
        future = TaskManager.taskPool.submit(task.run)
        TaskManager.futuresTaskMap[future] = task
        future.add_done_callback(taskCallback)
        #Configs.log("Future status {} for file {}".format(future._state, task.outputFile))
    else:
        TaskManager.observerTask = task
        TaskManager.observerSignal.set()

def taskCallback(future):
    task = TaskManager.futuresTaskMap.pop(future)
    task.future = future
    with TaskManager.managerLock:
        setTaskFinished(task)
        try:
            task.future.result()
            TaskManager.finishedTasks.add(task)
        except:
            TaskManager.failedTasks.add(task)
        TaskManager.runningTasks.remove(task)
        TaskManager.threadsUsed = TaskManager.threadsUsed - 1
        TaskManager.managerSignal.set()
        
def setTaskFinished(task):
    task.isFinished = True
    TaskManager.observerSignal.set()
        