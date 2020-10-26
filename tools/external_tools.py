'''
Created on Apr 14, 2020

@author: Vlad
'''

import subprocess
import os
import random
from configuration import Configs
from helpers import tasks

def generateMafftFilePathMap(inputPaths, outputDir):
    mafftMap = {inputPath : os.path.join(outputDir, "mafft_{}".format(os.path.basename(inputPath))) for inputPath in inputPaths}
    return mafftMap

def buildMafftAlignments(inputOutputPathMap):
    tasks = [buildMafftAlignment(inputPath, outputPath) for inputPath, outputPath in inputOutputPathMap.items()]
    return tasks
    
def buildMafftAlignment(inputPath, outputPath, subtablePath = None):
    return runMafft(inputPath, subtablePath, Configs.workingDir, outputPath, Configs.numCores)                

def runMafft(fastaPath, subtablePath, workingDir, outputPath, threads = 1):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    args = [Configs.mafftPath, "--localpair", "--maxiterate", "1000", "--ep", "0.123", 
            "--quiet", "--thread", str(threads), "--anysymbol"]
    if subtablePath is not None:
        args.extend(["--merge", subtablePath])
    args.extend([fastaPath, ">", tempPath])
    command = subprocess.list2cmdline(args)    
    copyMap = {tempPath : outputPath}
    return tasks.Task(command = command, workingDir = workingDir, outputFile = outputPath, fileCopyMap = copyMap)

def runMafftGuideTree(fastaPath, workingDir, outputPath, threads = 1):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    treeFile = os.path.join(os.path.dirname(fastaPath),  "{}.tree".format(os.path.basename(fastaPath)))
    args = [Configs.mafftPath, "--retree", "0", "--treeout", "--parttree",
            "--quiet", "--thread", str(threads), "--anysymbol"]
    args.extend(["--partsize", "1000"])
    args.extend([fastaPath, ">", tempPath])
    command = subprocess.list2cmdline(args)    
    copyMap = {treeFile : outputPath}
    return tasks.Task(command = command, workingDir = workingDir, outputFile = outputPath, fileCopyMap = copyMap)

def runMcl(matrixPath, inflation, workingDir, outputPath):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    args = [Configs.mclPath, matrixPath, "--abc", "-o", tempPath]
    if inflation is not None:
        args.extend(["-I", str(inflation)])
    command = subprocess.list2cmdline(args)
    copyMap = {tempPath : outputPath}
    return tasks.Task(command = command, workingDir = workingDir, outputFile = outputPath, fileCopyMap = copyMap)

def runMlrMcl(matrixPath, granularity, balance, inflation, workingDir, outputPath):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    args = [Configs.mlrmclPath, matrixPath, "-o", tempPath]
    if granularity is not None:
        args.extend(["-c", str(granularity)])
    if balance is not None:
        args.extend(["-b", str(balance)])    
    if inflation is not None:
        args.extend(["-i", str(inflation)])
    command = subprocess.list2cmdline(args)
    copyMap = {tempPath : outputPath}
    return tasks.Task(command = command, workingDir = workingDir, outputFile = outputPath, fileCopyMap = copyMap)

def runFastTree(fastaFilePath, workingDir, outputPath, mode = "normal", intree = None):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    
    args = [Configs.fasttreePath]
    if Configs.inferDataType(fastaFilePath) == "protein":
        args.extend(["-lg"])
    else:
        args.extend(["-nt", "-gtr"])
    
    if intree is not None:
        args.extend(["-intree", intree])
    
    if mode == "fast":
        #args.extend(["-fastest", "-quiet", "-mlnni", "4"]) 
        args.extend(["-fastest"]) 
    elif mode == "faster":
        #args.extend(["-fastest", "-quiet", "-mlnni", "4"]) 
        args.extend(["-fastest", "-mlnni", "4" ]) 
    elif mode == "fastest":
        #args.extend(["-fastest", "-quiet", "-noml"])
        args.extend(["-fastest", "-noml"])
    
    #logPath = os.path.join(workingDir, "fasttree_log.txt")
    #args.extend(["-log", logPath])    
    args.extend([fastaFilePath, ">", tempPath])
    command = subprocess.list2cmdline(args)
    copyMap = {tempPath : outputPath}
    return tasks.Task(command = command, workingDir = workingDir, outputFile = outputPath, fileCopyMap = copyMap)

def runRaxmlNg(fastaFilePath, workingDir, outputPath, threads = 8):
    # raxml-ng --msa prim.phy --model GTR+G --prefix T4 --threads 2 --seed 2 --tree pars{25},rand{25}
    baseName = os.path.basename(outputPath).replace(".","")
    raxmlFile = os.path.join(workingDir, "{}.raxml.bestTree".format(baseName))
    seed = random.randint(1, 1000000)
    args = [Configs.raxmlPath,
            "--msa", fastaFilePath,
            "--prefix", baseName,
            "--threads", str(threads),
            "--seed", str(seed)]
    
    if Configs.inferDataType(fastaFilePath) == "protein":
        args.extend(["--model", "LG+G"])
    else:
        args.extend(["--model", "GTR+G"])
        
    args.extend(["--tree", "pars{{{}}}".format(1)])
    
    command = subprocess.list2cmdline(args)
    copyMap = {raxmlFile : outputPath}
    return tasks.Task(command = command, workingDir = workingDir, outputFile = outputPath, fileCopyMap = copyMap)

def runHmmBuild(alignmentPath, workingDir, outputPath):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    args = [Configs.hmmbuildPath,'--ere', '0.59', "--cpu", "1"]
    args.extend(["--symfrac", "0.0", "--informat", "afa", tempPath, alignmentPath])
    command = subprocess.list2cmdline(args)
    copyMap = {tempPath : outputPath}
    return tasks.Task(command = command, workingDir = workingDir, outputFile = outputPath, fileCopyMap = copyMap)

def runHmmAlign(hmmModelPath, fragPath, workingDir, outputPath):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    args = [Configs.hmmalignPath, "-o", tempPath]
    args.extend([hmmModelPath, fragPath])
    command = subprocess.list2cmdline(args)
    copyMap = {tempPath : outputPath}
    return tasks.Task(command = command, workingDir = workingDir, outputFile = outputPath, fileCopyMap = copyMap)

def runHmmSearch(hmmModelPath, fragPath, workingDir, outputPath):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    args = [Configs.hmmsearchPath,"--noali", "--cpu", "1", "-o", tempPath, "-E", "99999999", "--max"]
    args.extend([hmmModelPath, fragPath])
    command = subprocess.list2cmdline(args)
    copyMap = {tempPath : outputPath}
    return tasks.Task(command = command, workingDir = workingDir, outputFile = outputPath, fileCopyMap = copyMap)