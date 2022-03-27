'''
Created on Oct 26, 2020

@author: Vlad
'''

import random
import threading
import os
import time
import concurrent.futures
from magus_configuration import Configs
from magus_tasks import files

'''
Launching and awaiting tasks.
To avoid deadlocks and stack overflows, only the main thread can submit tasks.
Thus, only the main thread runs alignment tasks, worker threads are used for other task types (like MAFFT).
'''

class TaskManager():
    
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
    serialTaskTypes = {"runAlignmentTask", "buildInducedSubalignment", "compressSubalignment"}
    contextStack = []
    
def startTaskManager():
    Configs.debug("Starting up the task manager..")
    
    tasksDir = os.path.join(Configs.workingDir, "tasks")
    TaskManager.pendingTasksDir = os.path.join(tasksDir, "tasks_pending")
    TaskManager.runningTasksFile = os.path.join(tasksDir, "tasks_running.txt")
    TaskManager.lockTasksFile = os.path.join(tasksDir, "tasks.lock")
    if not os.path.exists(TaskManager.pendingTasksDir):
        os.makedirs(TaskManager.pendingTasksDir)  
    
    TaskManager.managerPool = concurrent.futures.ThreadPoolExecutor(max_workers = 1)
    TaskManager.managerFuture = TaskManager.managerPool.submit(runTaskManager)
    TaskManager.taskPool = concurrent.futures.ThreadPoolExecutor(max_workers = Configs.numCores)
    Configs.debug("Task manager is up..")

def stopTaskManager():
    TaskManager.managerStopSignal = True
    with TaskManager.managerLock:
        TaskManager.managerSignal.set()
    try:
        Configs.debug("Winding down the task manager..")
        TaskManager.managerFuture.result()
    finally:
        Configs.log("Waiting for {} tasks to finish..".format(len(TaskManager.runningTasks)))
        TaskManager.taskPool.shutdown()
        dealWithFinishedTasks()        
        TaskManager.managerPool.shutdown()
        Configs.debug("Task manager stopped..")

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
                files.writeTasksToFile(stillRunningTasks, TaskManager.runningTasksFile, append = False)
        
    if len(TaskManager.failedTasks) > 0:
        processPendingTasks(TaskManager.failedTasks, 0, None)   
        
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
    
    if len(newTasks) > 0:
        launchedTasks, remainingTasks = processPendingTasks(newTasks, numToLaunch, None)
        numToLaunch = numToLaunch - len(launchedTasks)
    if numToLaunch > 0:
        pendingFiles = [os.path.join(TaskManager.pendingTasksDir, file) for file in os.listdir(TaskManager.pendingTasksDir) if file.endswith(".txt")]
        random.shuffle(pendingFiles)
        for taskFile in pendingFiles:
            if numToLaunch <= 0:
                break
            launchedTasks, remainingTasks = processPendingTasks(None, numToLaunch, taskFile)
            numToLaunch = numToLaunch - len(launchedTasks)

    TaskManager.submittedTasks = set()

def dealWithWaitingTasks():
    timeSinceFileCheck = time.time() - TaskManager.lastFilesCheckTime
    if timeSinceFileCheck >= 5: 
        for file, task in list(TaskManager.waitingTasks.items()):
            if os.path.exists(file):
                Configs.debug("Detected task completion: {}".format(file))
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

def processPendingTasks(tasks, numTasksToLaunch, taskFile):
    if taskFile is None:
        taskFile = os.path.join(TaskManager.pendingTasksDir, "task_file_{:03d}.txt".format(random.randint(0,999)))
    
    newTasks = True    
    with files.FileLock(taskFile.replace(".txt", ".lock")):
        if tasks is None:
            tasks = files.readTasksFromFile(taskFile)
            newTasks = False

        alignTasks = [t for t in tasks if t.taskType == "runAlignmentTask"]
        if len(alignTasks) > 0:
            with files.FileLock(TaskManager.lockTasksFile):
                runningTasks = set(files.readTasksFromFile(TaskManager.runningTasksFile))
                if len(runningTasks) > 0:
                    tasks = [t for t in tasks if t not in runningTasks]
                launchedTasks, remainingTasks = launchTasks(tasks, numTasksToLaunch)        
                writeRunningTasks = [t for t in launchedTasks if t.taskType == "runAlignmentTask"]
                files.writeTasksToFile(writeRunningTasks, TaskManager.runningTasksFile, append = True)
        else:
            launchedTasks, remainingTasks = launchTasks(tasks, numTasksToLaunch)
        
        files.writeTasksToFile(remainingTasks, taskFile, append = newTasks) 
        if len(launchedTasks) > 0:     
            Configs.debug("Launched {} tasks and deferred {} tasks..".format(len(launchedTasks), len(remainingTasks)))   
          
    return launchedTasks, remainingTasks      

def launchTasks(tasks, numTasksToLaunch):
    launchedTasks = []
    remainingTasks = []    
    threadsAvailable = Configs.numCores - TaskManager.threadsUsed
    toLaunch = min(numTasksToLaunch, threadsAvailable)    
    for task in tasks:
        if toLaunch > 0 and checkLaunchTask(task):
            Configs.debug("Launched a new task.. {}/{} threads used, type: {}, output file: {}".format(TaskManager.threadsUsed, Configs.numCores, task.taskType, task.outputFile))
            toLaunch = toLaunch - 1
            launchedTasks.append(task)
        else:
            remainingTasks.append(task)
    return launchedTasks, remainingTasks

def checkLaunchTask(task):
    if task.taskType not in TaskManager.serialTaskTypes:
        task.future = TaskManager.taskPool.submit(runTask, task)
        return True
    elif TaskManager.observerWaiting and TaskManager.observerTask is None:
        stack = TaskManager.contextStack
        if task.taskType != "runAlignmentTask" or len(stack) == 0 or task in stack[-1].subalignmentTasks:
            TaskManager.observerTask = task
            TaskManager.observerSignal.set()
            return True
    return False

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
                
        