'''
Created on Nov 1, 2020

@author: Vlad
'''

from tasks.manager import TaskManager, runTask
from configuration import Configs

'''
This is where the main thread goes to submit and await tasks.
When the main thread is in a blocking wait, it can pick up a new alignment task to run.
'''

def submitTasks(tasks):
    with TaskManager.managerLock:
        checkTaskManager()
        TaskManager.submittedTasks.update(tasks)
        TaskManager.managerSignal.set()

def asCompleted(tasks):
    unfinished = tasks
    while True:
        finished, unfinished = checkWhatFinished(unfinished)
        yield from finished
        if len(unfinished) == 0:
            return
        if len(finished) == 0:
            observeTaskManager()
        
def awaitTasks(tasks):
    finished, unfinished = checkWhatFinished(tasks)
    while len(unfinished) > 0:
        observeTaskManager()
        finished, unfinished = checkWhatFinished(unfinished)
        
def checkWhatFinished(tasks):
    finished, unfinished = [], []
    for t in tasks:
        finished.append(t) if t.checkFinished() else unfinished.append(t)
    return finished, unfinished

def observeTaskManager():
    TaskManager.observerWaiting = True
    TaskManager.managerSignal.set()
    TaskManager.observerSignal.wait(10)
    with TaskManager.managerLock:
        TaskManager.observerWaiting = False
        TaskManager.observerSignal.clear()
        checkTaskManager()
        task, TaskManager.observerTask = TaskManager.observerTask, None
    if task is not None:
        runTask(task)
            
        #manager.setTaskFinished(task)

def checkTaskManager():
    if not TaskManager.managerFuture.running():
        Configs.error("Task manager is dead for some reason..")
        TaskManager.managerFuture.result()
        raise Exception("Task manager is dead for some reason..")
        
