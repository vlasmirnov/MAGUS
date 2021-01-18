'''
Created on May 28, 2020

@author: Vlad
'''

import re
import os
import math
from tools import external_tools
from configuration import Configs
from helpers import sequenceutils

'''
def buildHmmScores(hmmPath, queriesPath, scoreMap):
    #tasks = [getHmmScores(hmmPath, queriesPath) for hmmPath in hmmPaths]
    queries = sequenceutils.readFromFasta(queriesPath, removeDashes = True)
    baseName = os.path.basename(queriesPath).split('.')[0]
    dirName = os.path.dirname(hmmPath)
    
    taxonScores = {}
    for taxon in queries:
        inputName = os.path.join(dirName, "{}.txt".format(baseName))
        outputName = os.path.join(dirName, "{}_score.txt".format(baseName))
        if os.path.exists(outputName):
            os.remove(outputName)
        sequenceutils.writeFasta({taxon : queries[taxon]}, inputName)
        getHmmScores(hmmPath, inputName, outputName).run()
        subsetScores = readSearchFile(outputName)
        taxonScores[taxon] = subsetScores[taxon][1]
    scoreMap[hmmPath] = taxonScores
'''

def buildHmmScores(hmmPaths, queriesPath, scoreFileHmmFileMap):
    #tasks = [getHmmScores(hmmPath, queriesPath) for hmmPath in hmmPaths]
    queries = sequenceutils.readFromFasta(queriesPath, removeDashes = True)
    baseName = os.path.basename(queriesPath).split('.')[0]
    dirName = os.path.join(os.path.dirname(queriesPath), "chunks_{}".format(baseName))
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    
    chunkSize = 1000
    
    taxa = list(queries.keys())
    inputOutputs = []
    for i in range(math.ceil(len(taxa) / chunkSize)):
        chunk = taxa[i*chunkSize : min(len(taxa), (i+1)*chunkSize)]
        inputName = os.path.join(dirName, "{}_chunk_{}.txt".format(baseName, i+1))
        sequenceutils.writeFasta(queries, inputName, chunk)
        for hmmPath in hmmPaths:
            outputName = os.path.join(os.path.dirname(hmmPath), "{}_chunk_{}_score.txt".format(baseName, i+1))
            inputOutputs.append((hmmPath, inputName, outputName))
            scoreFileHmmFileMap[outputName] = hmmPath
    
    tasks = [getHmmScores(hmmPath, inputPath, outputPath) for hmmPath, inputPath, outputPath in inputOutputs]
    return tasks

def getHmmScores(hmmPath, queriesPath, scorePath):
    workingDir = os.path.dirname(hmmPath)
    #searchPath = os.path.join(workingDir, "hmm_search.txt")
    task = external_tools.runHmmSearch(hmmPath, queriesPath, workingDir, scorePath)
    return task

def readHmmScores(searchFiles):
    sequenceScores = {}
    for file in searchFiles:
        subsetScores = readSearchFile(file)
        taxonScores = {taxon : scores[1] for taxon, scores in subsetScores.items()}
        sequenceScores[file] = taxonScores
    return sequenceScores
    
def buildHmms(sequencesHmmsPathsMap):
    tasks = [buildHmmOverAlignment(sequencePath, hmmPath) for sequencePath, hmmPath in sequencesHmmsPathsMap.items()]
    return tasks
    
def buildHmmOverAlignment(sequencePath, hmmPath):
    workingDir = os.path.dirname(hmmPath)
    task = external_tools.runHmmBuild(sequencePath, workingDir, hmmPath)
    return task

def combineHmmAlignments(alignFiles, outputAlignmentPath, includeInsertions):
    alignment = {}
    for file in alignFiles:
        alignment.update(sequenceutils.readFromStockholm(file, includeInsertions))
    sequenceutils.writeFasta(alignment, outputAlignmentPath, None)

def mergeHmmAlignments(alignFiles, outputAlignmentPath, includeInsertions):
    for file in alignFiles:
        alignment = sequenceutils.readFromStockholm(file, includeInsertions)
        sequenceutils.writeFasta(alignment, outputAlignmentPath, None, True)

def hmmAlignQueries(hmmPath, queriesPath):
    queries = sequenceutils.readFromFasta(queriesPath, removeDashes = True)
    baseName = os.path.basename(queriesPath).split('.')[0]
    dirName = os.path.join(os.path.dirname(queriesPath), "chunks_{}".format(baseName))
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    chunkSize = 1000
    
    taxa = list(queries.keys())
    alignFiles = {}
    for i in range(math.ceil(len(taxa) / chunkSize)):
        chunk = taxa[i*chunkSize : min(len(taxa), (i+1)*chunkSize)]
        inputName = os.path.join(dirName, "{}_chunk_{}.txt".format(baseName, i+1))
        outputName = os.path.join(dirName, "{}_chunk_{}_aligned.txt".format(baseName, i+1))
        sequenceutils.writeFasta(queries, inputName, chunk)
        alignFiles[inputName] = outputName
    
    tasks = []
    for inputPath, outputPath in alignFiles.items():
        task = buildHmmAlignment(hmmPath, inputPath, outputPath)
        tasks.append(task)
    return tasks

def buildHmmAlignment(hmmPath, queriesPath, outputAlignmentPath):
    workingDir = os.path.dirname(hmmPath)
    task = external_tools.runHmmAlign(hmmPath, queriesPath, workingDir, outputAlignmentPath)
    return task

#from PASTA repo
def readSearchFile(searchFilePath):
    with open(searchFilePath, 'r') as searchFile:
        results = {}

        pattern = re.compile(
            r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+"
            r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)")
        start_reading = False
        for line in searchFile:
            line = line.strip()
            if (not start_reading and line.startswith("E-value") is True):
                start_reading = True
            elif (start_reading and line == ""):
                start_reading = False
                break
            elif (start_reading):
                matches = pattern.search(line)
                if (matches is not None and matches.group(0).find("--") == -1):
                    results[matches.group(9).strip()] = (
                        float(matches.group(1).strip()),
                        float(matches.group(2).strip()))
                    # _LOG.debug("Fragment scores;"
                    #           "fragment:%s E-Value:%s BitScore:%s" %(matches
                    # .group(9).strip(),matches.group(1).strip(), matches.
                    # group(2).strip()))
        return results