'''
Created on May 28, 2020

@author: Vlad
'''

import re
import os
import shutil
import concurrent.futures
from tools import external_tools

def buildHmmScores(sequencesHmmsPathsMap, queriesPath):
    print("Launching {} HMM scoring jobs with {} workers..".format(len(sequencesHmmsPathsMap), os.cpu_count()))
    sequenceScores = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers = os.cpu_count()) as executor:
        jobs = {executor.submit(getHmmScores, sequencesHmmsPathsMap[sequencePath], queriesPath) : sequencePath
                   for sequencePath in sequencesHmmsPathsMap}
        for job in concurrent.futures.as_completed(jobs):
            try:
                taxonScores = job.result()
                sequenceScores[jobs[job]] = taxonScores
            except Exception as exc:
                print("Worker for HMM scores {} threw an exception:\n{}".format(jobs[job], exc))
                raise
    print("All HMM scores workers done..")
    return sequenceScores

def getHmmScores(hmmDir, queriesPath):
    searchPath = os.path.join(hmmDir, "hmm_search.txt")
    hmmPath = os.path.join(hmmDir, "hmm_model.txt")
    external_tools.runHmmSearch(hmmPath, queriesPath, hmmDir, searchPath)
    subsetScores = readSearchFile(searchPath)
    taxonScores = {}
    for taxon, scores in subsetScores.items():
        taxonScores[taxon] = scores[1]
    return taxonScores

def buildHmms(sequencesHmmsPathsMap):
    print("Launching {} HMM builds with {} workers..".format(len(sequencesHmmsPathsMap), os.cpu_count()))
    
    with concurrent.futures.ThreadPoolExecutor(max_workers = os.cpu_count()) as executor:
        jobs = {executor.submit(buildHmmOverAlignment, sequencesHmmsPathsMap[sequencePath], sequencePath) : sequencePath
                   for sequencePath in sequencesHmmsPathsMap}
        for job in concurrent.futures.as_completed(jobs):
            try:
                job.result()
            except Exception as exc:
                print("Worker for HMM {} threw an exception:\n{}".format(jobs[job], exc))
                raise
    print("All HMM workers done..")
    
def buildHmmOverAlignment(hmmDir, alignmentPath):
    if os.path.exists(hmmDir):
        shutil.rmtree(hmmDir)
    os.makedirs(hmmDir)  
    hmmPath = os.path.join(hmmDir, "hmm_model.txt")
    external_tools.runHmmBuild(alignmentPath, hmmDir, hmmPath)
    return hmmPath

def hmmAlignQueries(hmmDir, queriesPath):
    hmmPath = os.path.join(hmmDir, "hmm_model.txt")
    alignResultPath = os.path.join(hmmDir, "hmm_align_result.txt")
    external_tools.runHmmAlign(hmmPath, queriesPath, hmmDir, alignResultPath)
    return alignResultPath

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