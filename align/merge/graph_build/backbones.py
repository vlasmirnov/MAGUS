'''
Created on Jul 28, 2021

@author: Vlad
'''

import os
import random
import heapq
import math

from helpers import sequenceutils
from tasks import task
from configuration import Configs
from tools import external_tools


def requestMafftBackbones(context):
    missingBackboneFiles = {}    
    for n in range(Configs.mafftRuns):
        unalignedFile = os.path.join(context.graph.workingDir, "backbone_{}_unalign.txt".format(n+1))
        alignedFile = os.path.join(context.graph.workingDir, "backbone_{}_mafft.txt".format(n+1))
        if os.path.exists(alignedFile):
            Configs.log("Existing backbone file found: {}".format(alignedFile))            
            context.backbonePaths.append(alignedFile)
        else:
            missingBackboneFiles[unalignedFile] = alignedFile
            
        if Configs.graphBuildHmmExtend:
            context.backboneExtend.add(alignedFile)

    assignBackboneTaxa(context, missingBackboneFiles)
    for unalignedFile, alignedFile in missingBackboneFiles.items():
        backboneTask = external_tools.buildMafftAlignment(unalignedFile, alignedFile)
        context.backboneTasks.append(backboneTask)
    
    if not Configs.graphBuildHmmExtend:
        for file in list(missingBackboneFiles.keys()) + context.backbonePaths:
            backbone = sequenceutils.readFromFasta(file)
            context.backboneTaxa.update(backbone)        
            
    task.submitTasks(context.backboneTasks)    
    
def assignBackboneTaxa(context, missingBackbones):
    if len(missingBackbones) == 0:
        return
    
    numTaxa = max(1, int(Configs.mafftSize/len(context.subsetPaths)))
    backbones = {file : {} for file in missingBackbones}
    
    if Configs.graphBuildStrategy.lower() == "random":    
        buildBackbonesRandom(context, backbones, numTaxa)
        
    elif Configs.graphBuildStrategy.lower() == "longest":
        buildBackbonesLongest(context, backbones, numTaxa)
    
    elif Configs.graphBuildStrategy.lower() == "longestrandom":
        buildBackbonesLongestRandom(context, backbones, numTaxa)
    
    elif Configs.graphBuildStrategy.lower() == "longestrandom2":
        buildBackbonesLongestRandom2(context, backbones, numTaxa)
    
    elif Configs.graphBuildStrategy.lower() == "longestrandom3":
        buildBackbonesLongestRandom3(context, backbones, numTaxa)
    
    elif Configs.graphBuildStrategy.lower() == "longestrandom4":
        buildBackbonesLongestRandom4(context, backbones, numTaxa)
    
    elif Configs.graphBuildStrategy.lower() == "longestrandom5":
        buildBackbonesLongestRandom5(context, backbones, numTaxa)
    
    elif Configs.graphBuildStrategy.lower() == "longestrandom6":
        buildBackbonesLongestRandom6(context, backbones, numTaxa)
    
    elif Configs.graphBuildStrategy.lower() == "longestrandom7":
        buildBackbonesLongestRandom7(context, backbones, numTaxa)

    elif Configs.graphBuildStrategy.lower() == "coverage":
        buildBackbonesCoverage(context, backbones, numTaxa)  
    
    elif Configs.graphBuildStrategy.lower() == "coverage2":
        buildBackbonesCoverage2(context, backbones, numTaxa)                          
                    
    for file, backbone in backbones.items():
        sequenceutils.writeFasta(backbone, file)
    
    
def buildBackbonesRandom(context, backbones, numTaxa):
    Configs.log("Preparing {} backbones with {} RANDOM sequences per subset..".format(len(backbones), numTaxa))
    for file, backbone in backbones.items():
        for subset in context.subsets:
            random.shuffle(subset)
            for taxon in subset[:numTaxa]:
                backbone[taxon] = context.unalignedSequences[taxon]       

def buildBackbonesLongest(context, backbones, numTaxa):
    Configs.log("Preparing {} backbones with {} LONGEST sequences per subset..".format(len(backbones), numTaxa))
    for subset in context.subsets:
        sortedByLength = sorted(subset, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        curIdx = 0
        for file, backbone in backbones.items():
            window = sortedByLength[curIdx : curIdx + numTaxa]
            if curIdx + numTaxa > len(sortedByLength):
                window.extend(sortedByLength[0 : min(curIdx, curIdx + numTaxa - len(sortedByLength))])
            for taxon in window:
                backbone[taxon] = context.unalignedSequences[taxon]
            curIdx = (curIdx + numTaxa) % len(sortedByLength)

def buildBackbonesLongestRandom(context, backbones, numTaxa):
    Configs.log("Preparing {} backbones with {} LONGEST (randomized) sequences per subset..".format(len(backbones), numTaxa))
    for subset in context.subsets:
        sortedByLength = sorted(subset, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        topSequences = sortedByLength[ : max(numTaxa, int(0.25*len(sortedByLength)))]
        for file, backbone in backbones.items():
            random.shuffle(topSequences)
            for taxon in topSequences[:numTaxa]:
                backbone[taxon] = context.unalignedSequences[taxon]

def buildBackbonesLongestRandom2(context, backbones, numTaxa):
    Configs.log("Preparing {} backbones with {} LONGEST (randomized-2) sequences per subset..".format(len(backbones), numTaxa))
    for subset in context.subsets:
        sortedByLength = sorted(subset, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        topSequences = sortedByLength[ : numTaxa * 2]
        for file, backbone in backbones.items():
            random.shuffle(topSequences)
            for taxon in topSequences[:numTaxa]:
                backbone[taxon] = context.unalignedSequences[taxon]

def buildBackbonesLongestRandom3(context, backbones, numTaxa):
    Configs.log("Preparing {} backbones with {} LONGEST (randomized-3) sequences per subset..".format(len(backbones), numTaxa))
    for subset in context.subsets:
        sortedByLength = sorted(subset, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        topSequences = sortedByLength[ : numTaxa * 3]
        for file, backbone in backbones.items():
            random.shuffle(topSequences)
            for taxon in topSequences[:numTaxa]:
                backbone[taxon] = context.unalignedSequences[taxon]

def buildBackbonesLongestRandom4(context, backbones, numTaxa):
    Configs.log("Preparing {} backbones with {} LONGEST (randomized-4) sequences per subset..".format(len(backbones), numTaxa))
    for subset in context.subsets:
        sortedByLength = sorted(subset, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        numTop = math.ceil(numTaxa/2)
        topSeq, remainder = sortedByLength[ : numTop], sortedByLength[numTop : ]
        for file, backbone in backbones.items():
            random.shuffle(remainder)
            for taxon in topSeq + remainder[ : numTaxa-numTop]:
                backbone[taxon] = context.unalignedSequences[taxon]

def buildBackbonesLongestRandom5(context, backbones, numTaxa):
    Configs.log("Preparing {} backbones with {} LONGEST (randomized-5) sequences per subset..".format(len(backbones), numTaxa))
    for subset in context.subsets:
        sortedByLength = sorted(subset, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        numTop = math.ceil(numTaxa/4)
        topSeq, remainder = sortedByLength[ : numTop], sortedByLength[numTop : ]
        for file, backbone in backbones.items():
            random.shuffle(remainder)
            for taxon in topSeq + remainder[ : numTaxa-numTop]:
                backbone[taxon] = context.unalignedSequences[taxon]

def buildBackbonesLongestRandom6(context, backbones, numTaxa):
    Configs.log("Preparing {} backbones with {} LONGEST (randomized-6) sequences per subset..".format(len(backbones), numTaxa))
    for subset in context.subsets:
        sortedByLength = sorted(subset, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        #topQuartileTaxon = sortedByLength[int(0.25*len(sortedByLength))]
        topLength = len(context.unalignedSequences[sortedByLength[0]].seq)
        fullLength = []
        notFullLength = []
        for t in sortedByLength:
            if abs(len(context.unalignedSequences[t].seq) - topLength) < 0.25 * topLength:
                fullLength.append(t)
            else:
                notFullLength.append(t) 
        
        Configs.log("Found {}/{} full-length sequences..".format(len(fullLength), len(subset)))        
        for file, backbone in backbones.items():
            random.shuffle(fullLength)
            random.shuffle(notFullLength)  
            allTaxa = fullLength + notFullLength  
            for taxon in allTaxa[:numTaxa]:
                backbone[taxon] = context.unalignedSequences[taxon]

def buildBackbonesLongestRandom7(context, backbones, numTaxa):
    Configs.log("Preparing {} backbones with {} LONGEST (randomized-7) sequences per subset..".format(len(backbones), numTaxa))
    for subset in context.subsets:
        sortedByLength = sorted(subset, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        #topQuartileTaxon = sortedByLength[int(0.25*len(sortedByLength))]
        topLength = len(context.unalignedSequences[sortedByLength[0]].seq)
        fullLength = []
        notFullLength = []
        for t in sortedByLength:
            if abs(len(context.unalignedSequences[t].seq) - topLength) < 0.5 * topLength:
                fullLength.append(t)
            else:
                notFullLength.append(t) 
        
        Configs.log("Found {}/{} full-length sequences..".format(len(fullLength), len(subset)))        
        for file, backbone in backbones.items():
            random.shuffle(fullLength)
            random.shuffle(notFullLength)  
            allTaxa = fullLength + notFullLength  
            for taxon in allTaxa[:numTaxa]:
                backbone[taxon] = context.unalignedSequences[taxon]

def buildBackbonesCoverage(context, backbones, numTaxa):
    context.awaitSubalignments()
    Configs.log("Preparing {} backbones with {} COVERAGE sequences per subset..".format(len(backbones), numTaxa))
        
    for subalignPath in context.subalignmentPaths:
        subalignment = sequenceutils.readFromFasta(subalignPath, removeDashes=False)
        taxa = list(subalignment.keys())
        sortedByLength = sorted(taxa, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        topSequences = sortedByLength[ : numTaxa*2]
        
        coverArray = None
        coverageDict = {}
        for taxon in sortedByLength:
            sequence = subalignment[taxon]
            if coverArray is None:
                coverArray = [[] for i in range(len(sequence.seq))]
            coverageDict[taxon] = []
            for i, c in enumerate(sequence.seq):
                if c not in ('-', '.', '_'):
                    coverArray[i].append(taxon)
                    coverageDict[taxon].append(i)
        
        
        for file, backbone in backbones.items():
            heap = []
            usedTaxons = set()
            coverage = [0] * len(coverArray)
            coverArrayPointers = [0] * len(coverArray)
            
            random.shuffle(topSequences)
            numSeq = int(numTaxa / 2)
            for taxon in topSequences[:numSeq]:
                usedTaxons.add(taxon)
                backbone[taxon] = context.unalignedSequences[taxon]
                for pos in coverageDict[taxon]:
                    coverage[pos] = coverage[pos] + 1
                
                
            for i in range(len(coverArray)):
                for p in range(len(coverArray[i])):
                    taxon = coverArray[i][p]
                    if not taxon in usedTaxons:
                        coverArrayPointers[i] = p
                        item = (coverage[i], -1 * len(coverageDict[taxon]), random.random(), i, taxon)
                        heapq.heappush(heap, item)
                        break
            
            while len(usedTaxons) < min(numTaxa, len(subalignment)):
                colCover, taxCover, trash, i, taxon = heapq.heappop(heap)
                if taxon in usedTaxons:
                    continue
                if colCover < coverage[i]:
                    item = (coverage[i], -1 * len(coverageDict[taxon]), random.random(), i, taxon)
                    heapq.heappush(heap, item)
                    continue
                
                usedTaxons.add(taxon)
                backbone[taxon] = context.unalignedSequences[taxon]                
                for pos in coverageDict[taxon]:
                    coverage[pos] = coverage[pos] + 1
                    
                    if taxon == coverArray[pos][coverArrayPointers[pos]]:
                        while coverArrayPointers[pos] < len(coverArray[pos]) and coverArray[pos][coverArrayPointers[pos]] in usedTaxons:
                            coverArrayPointers[pos] = coverArrayPointers[pos] + 1
                        if coverArrayPointers[pos] < len(coverArray[pos]):
                            newTaxon = coverArray[pos][coverArrayPointers[pos]]
                            item = (coverage[pos], -1 * len(coverageDict[newTaxon]), random.random(), pos, newTaxon)
                            heapq.heappush(heap, item) 

def buildBackbonesCoverage2(context, backbones, numTaxa):
    context.awaitSubalignments()
    Configs.log("Preparing {} backbones with {} COVERAGE2 sequences per subset..".format(len(backbones), numTaxa))
        
    for subalignPath in context.subalignmentPaths:
        subalignment = sequenceutils.readFromFasta(subalignPath, removeDashes=False)
        taxa = list(subalignment.keys())
        sortedByLength = sorted(taxa, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        
        coverArray = None
        coverageDict = {}
        for taxon in sortedByLength:
            sequence = subalignment[taxon]
            if coverArray is None:
                coverArray = [[] for i in range(len(sequence.seq))]
            coverageDict[taxon] = []
            for i, c in enumerate(sequence.seq):
                if c not in ('-', '.', '_'):
                    coverArray[i].append(taxon)
                    coverageDict[taxon].append(i)
        
        
        for file, backbone in backbones.items():
            heap = []
            heapTaxons = set()
            usedTaxons = set()
            usedSites = set()
             
            for i in range(len(coverArray)):
                for p in range(len(coverArray[i])):
                    taxon = coverArray[i][p]
                    if taxon not in heapTaxons:
                        item = (-1 * len(coverageDict[taxon]), random.random(), taxon)
                        heapTaxons.add(taxon)
                        heapq.heappush(heap, item)
                        break
            
            while len(usedTaxons) < min(numTaxa, len(subalignment)) and len(heap) > 0:
                taxCover, trash, taxon = heapq.heappop(heap)   
                for pos in coverageDict[taxon]:
                    if pos not in usedSites:
                        usedSites.add(pos)
                        usedTaxons.add(taxon)
                        backbone[taxon] = context.unalignedSequences[taxon]                
            Configs.log("Selected {} coverage sequences..".format(len(usedTaxons)))
            
            remainder = [t for t in taxa if t not in usedTaxons]
            random.shuffle(remainder)
            for taxon in remainder[ : numTaxa - len(usedTaxons)]:
                backbone[taxon] = context.unalignedSequences[taxon]