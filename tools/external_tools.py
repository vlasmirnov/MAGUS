'''
Created on Apr 14, 2020

@author: Vlad
'''

import subprocess
import os
import shutil
import random
import concurrent.futures
from configuration import Configs


def generateMafftFilePathMap(inputPaths, outputDir):
    mafftMap = {}
    for inputPath in inputPaths:
        alignedPath = os.path.join(outputDir, "mafft_{}".format(os.path.basename(inputPath)))
        mafftMap[inputPath] = alignedPath
    return mafftMap

def buildMafftAlignments(inputOutputPathMap):
    print("Launching {} MAFFT alignments with {} workers..".format(len(inputOutputPathMap), Configs.numCores))
    
    with concurrent.futures.ThreadPoolExecutor(max_workers = Configs.numCores) as executor:
        jobs = {executor.submit(buildMafftAlignment, inputPath, inputOutputPathMap[inputPath]) : inputPath
                   for inputPath in inputOutputPathMap}
        for job in concurrent.futures.as_completed(jobs):
            try:
                job.result()
            except Exception as exc:
                print("Worker for MAFFT alignment {} threw an exception:\n{}".format(jobs[job], exc))
                raise
    print("All MAFFT alignment workers done..")

def buildMafftAlignment(inputPath, outputPath, subtablePath = None):
    if not os.path.exists(outputPath):         
        print("Launching MAFFT on {}..".format(inputPath))  
        runMafft(inputPath, subtablePath, Configs.workingDir, outputPath, Configs.numCores)                
    print("Completed MAFFT on {}..".format(inputPath))

def runMafft(fastaPath, subtablePath, workingDir, outputPath, threads = 1):
    args = [Configs.mafftPath, "--localpair", "--maxiterate", "1000", "--ep", "0.123", 
            "--quiet", "--thread", str(threads), "--anysymbol"]
    if subtablePath is not None:
        args.extend(["--merge", subtablePath])
    args.extend([fastaPath, ">", outputPath])
    command = subprocess.list2cmdline(args)
    Configs.log("Running command: {}".format(command))
    subprocess.run(command, shell=True, cwd = workingDir)
    return command

def runMcl(matrixPath, inflation, workingDir, outputPath):
    args = [Configs.mclPath, matrixPath, "--abc", "-o", outputPath]
    if inflation is not None:
        args.extend(["-I", str(inflation)])
    command = subprocess.list2cmdline(args)
    Configs.log("Running command: {}".format(command))
    subprocess.run(command, shell=True, cwd = workingDir)
    return command

def runFastTree(fastaFilePath, workingDir, outputPath):
    #tempFile = os.path.join(workingDir, os.path.basename(outputPath))
    args = [Configs.fasttreePath]
    if Configs.inferDataType(fastaFilePath) == "protein":
        args.extend(["-lg"])
    else:
        args.extend(["-nt", "-gtr"])
    args.extend([fastaFilePath, ">", outputPath])
    command = subprocess.list2cmdline(args)
    Configs.log("Running command: {}".format(command))
    subprocess.run(command, shell=True, cwd = workingDir)
    return command

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
    Configs.log("Running command: {}".format(command))
    subprocess.run(command, shell=True, cwd = workingDir)
    shutil.copyfile(raxmlFile, outputPath)
    return command

def runHmmBuild(alignmentPath, workingDir, outputPath):
    #tempFile = os.path.join(workingDir, os.path.basename(outputPath))    
    #args = [hmmBuildPath,'--ere', '0.59', "--cpu", "1", "--dna"]
    args = [Configs.hmmbuildPath,'--ere', '0.59', "--cpu", "1"]
    args.extend(["--symfrac", "0.0", "--informat", "afa", outputPath, alignmentPath])
    command = subprocess.list2cmdline(args)
    Configs.log("Running command: {}".format(command))
    subprocess.run(command, shell=True, cwd = workingDir)
    return command

def runHmmAlign(hmmModelPath, fragPath, workingDir, outputPath):
    #tempFile = os.path.join(workingDir, os.path.basename(outputPath))    
    #args = [hmmAlignPath,"--dna", "-o", outputPath]
    args = [Configs.hmmalignPath, "-o", outputPath]
    args.extend([hmmModelPath, fragPath])
    command = subprocess.list2cmdline(args)
    Configs.log("Running command: {}".format(command))
    subprocess.run(command, shell=True, cwd = workingDir)
    return command

def runHmmSearch(hmmModelPath, fragPath, workingDir, outputPath):
    #tempFile = os.path.join(workingDir, os.path.basename(outputPath))    
    args = [Configs.hmmsearchPath,"--noali", "--cpu", "1", "-o", outputPath, "-E", "99999999", "--max"]
    args.extend([hmmModelPath, fragPath])
    command = subprocess.list2cmdline(args)
    Configs.log("Running command: {}".format(command))
    subprocess.run(command, shell=True, cwd = workingDir)
    return command