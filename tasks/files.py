'''
Created on Oct 26, 2020

@author: Vlad
'''

import os
import json
import time
import random
from tasks import task


def writeTasksToFile(taskList, tasksFile, append = True):
    if not append and len(taskList) == 0 and os.path.exists(tasksFile):
        os.remove(tasksFile)
    elif len(taskList) > 0:
        with open(tasksFile, 'a' if append else 'w') as file:
            for t in taskList:
                file.write(t.json + "\n")

def readTasksFromFile(tasksFile):
    fileTasks = []
    if os.path.exists(tasksFile):
        with open(tasksFile) as file:
            for line in file:
                mapper = json.loads(line.strip())
                fileTasks.append(task.Task(**mapper))
    return fileTasks


class FileLock:
    
    def __init__(self, filePath):
        self.filePath = filePath
    
    def __enter__(self):
        while True:
            try:
                lock = open(self.filePath, 'x')
                lock.close()
                return self
            except:
                time.sleep(random.random()*0.1 + 0.05)
                #time.sleep(random.random() + 0.5)
            
    def __exit__(self, excType, excVal, excTb):
        os.remove(self.filePath)
