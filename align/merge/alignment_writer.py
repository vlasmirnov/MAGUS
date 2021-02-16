'''
Created on Jan 6, 2021

@author: Vlad
'''

import os
import shutil
import heapq
from configuration import Configs
from helpers import sequenceutils
from tasks import task


def writeAlignment(context):
    context.graph.addSingletonClusters()
    if Configs.constrain:
        compressAlignment(context)
        writeUnpackedAlignment(context)
    else:
        writeUnconstrainedAlignment(context) 

def compressAlignment(context):
    graph = context.graph
    numCells = len(context.unalignedSequences) * len(graph.clusters)
    Configs.log("Uncompressed alignment will have {} cells..".format(numCells))
    if numCells > Configs.alignmentSizeLimit*1e9:
        Configs.log("Alignment will be more than {} Gigs, compressing..".format(Configs.alignmentSizeLimit))
        compressions, numLetters = buildCompressions(context)
        startCols, startHoms = len(graph.clusters), countHomologies(graph.clusters, numLetters, set())
                
        compressedClusters, insertions = compressClusters(context, graph.clusters, compressions, numLetters, Configs.alignmentSizeLimit*1e9)
        
        endCols, endHoms = len(compressedClusters), countHomologies(compressedClusters, numLetters, insertions)
        Configs.log("Uncompressed alignment has {} columns and {} homologous pairs..".format(startCols, startHoms))
        Configs.log("Compressed alignment has {} columns and {} homologous pairs..".format(endCols, endHoms))
        Configs.log("{}% columns and {}% homologous pairs lost..".format(100*(startCols - endCols)/startCols, 100*(startHoms - endHoms)/startHoms))
        numCells = len(context.unalignedSequences) * len(compressedClusters)
        Configs.log("Compressed alignment will have {} cells..".format(numCells))
        
        graph.clusters = compressedClusters
        graph.insertions = insertions 
    
def buildInducedSubalignment(**kwargs):
    alignmentColumnsPath = kwargs["alignmentColumnsPath"]
    subalignmentPath = kwargs["subalignmentPath"]
    inducedAlignPath = kwargs["outputFile"]
    tempInducedAlignPath = os.path.join(os.path.dirname(inducedAlignPath), "temp_{}".format(os.path.basename(inducedAlignPath)))
    
    alignColumns = []
    with open(alignmentColumnsPath) as f:
        insertIdxs = set([int(token) for token in f.readline().strip().split()])
        for line in f:
            tokens = set([int(token) for token in line.strip().split()])
            alignColumns.append(tokens) 
    
    subsetAlign = sequenceutils.readFromFasta(subalignmentPath, removeDashes=False)
    inducedAlign = {taxon : ['-'] * len(alignColumns) for taxon in subsetAlign}
        
    for idx, column in enumerate(alignColumns):
        for taxon in subsetAlign:
            for c in column:
                letter = subsetAlign[taxon].seq[c]
                if letter != '-':
                    letter = letter.lower() if c in insertIdxs else letter
                    inducedAlign[taxon][idx] = letter

    for s in inducedAlign:
        inducedAlign[s] = sequenceutils.Sequence(s, "".join(inducedAlign[s]))
    sequenceutils.writeFasta(inducedAlign, tempInducedAlignPath)
    shutil.move(tempInducedAlignPath, inducedAlignPath)
    
def writeUnpackedAlignment(context):
    graph = context.graph
    filePath = context.outputFile
    
    tempFile = os.path.join(os.path.dirname(filePath), "temp_{}".format(os.path.basename(filePath)))
    if os.path.exists(tempFile):
        os.remove(tempFile)
        
    clusterMap = {path : [[] for c in graph.clusters] for path in context.subalignmentPaths}
    for idx, cluster in enumerate(graph.clusters):
        for b in cluster:
            bsub, bpos = graph.matSubPosMap[b] 
            clusterMap[context.subalignmentPaths[bsub]][idx].append(bpos)
    
    inserts = {path : [] for path in context.subalignmentPaths}
    for b in graph.insertions:
        bsub, bpos = graph.matSubPosMap[b] 
        inserts[context.subalignmentPaths[bsub]].append(bpos)
        
    Configs.log("Assembling final alignment in {}".format(filePath))
    inducedSubalignTasks = []
    for subalignPath in context.subalignmentPaths:
        alignmentColumnsPath = os.path.join(context.graph.workingDir, "alignment_columns_{}".format(os.path.basename(subalignPath)))
        with open(alignmentColumnsPath, 'w') as textFile:
            textFile.write("{}\n".format(" ".join([str(c) for c in inserts[subalignPath]])))
            for cluster in clusterMap[subalignPath]:
                textFile.write("{}\n".format(" ".join([str(c) for c in cluster])))
        
        inducedAlignPath = os.path.join(graph.workingDir, "induced_{}".format(os.path.basename(subalignPath)))
        args = {"alignmentColumnsPath" : alignmentColumnsPath, "subalignmentPath" : subalignPath, "outputFile" : inducedAlignPath}
        inducedTask = task.Task(taskType = "buildInducedSubalignment", outputFile = args["outputFile"], taskArgs = args)
        inducedSubalignTasks.append(inducedTask)
        #inducedTask.submitTask()
    
    task.submitTasks(inducedSubalignTasks)            
    for inducedTask in task.asCompleted(inducedSubalignTasks):
        inducedAlign = sequenceutils.readFromFasta(inducedTask.outputFile, removeDashes=False)
        Configs.log("Appending induced alignment, {} sequences of length {}..".format(len(inducedAlign), len(next(iter(inducedAlign.values())).seq)))
        sequenceutils.writeFasta(inducedAlign, tempFile, append = True)   
        
        os.remove(inducedTask.taskArgs["alignmentColumnsPath"])
        os.remove(inducedTask.outputFile)
    shutil.move(tempFile, filePath)
    Configs.log("Wrote final alignment to {}".format(filePath))        
    #Configs.log("Wrote out {} clusters..".format(len(graph.clusters)))

def buildCompressions(context):
    graph = context.graph
    idxMap = {}
    for cluster in graph.clusters:
        for b in cluster:
            bsub, bpos = graph.matSubPosMap[b] 
            idxMap[context.subalignmentPaths[bsub], bpos] = b
    
    compressionTasks = []    
    for subalignPath in context.subalignmentPaths:        
        compressPath = os.path.join(context.graph.workingDir, "compression_{}".format(os.path.basename(subalignPath)))
        args = {"subalignmentPath" : subalignPath, "outputFile" : compressPath}
        compressionTask = task.Task(taskType = "compressSubalignment", outputFile = args["outputFile"], taskArgs = args)
        compressionTasks.append(compressionTask)
    task.submitTasks(compressionTasks)
    
    compressions = {}
    numLetters = {}
    for compressionTask in task.asCompleted(compressionTasks):
        with open(compressionTask.outputFile) as f:
            letterCounts = [int(token) for token in f.readline().strip().split()]
            for i, letters in enumerate(letterCounts):
                numLetters[idxMap[compressionTask.taskArgs["subalignmentPath"], i]] = letters
            comps = [int(token) for token in f.readline().strip().split()]
            for i, dest in enumerate(comps):
                if dest != -1:
                    compressions[idxMap[compressionTask.taskArgs["subalignmentPath"], i]] = idxMap[compressionTask.taskArgs["subalignmentPath"], dest]
        #os.remove(compressionTask.taskArgs["singletonsPath"])
        #os.remove(compressionTask.outputFile)     
    return compressions, numLetters

def combineCompressions(context, compressions):
    graph = context.graph
    
    mergedClusters = []
    compMap = {}
    for cluster in graph.clusters:
        newCluster = []
        for b in cluster:
            if b in compressions:
                center = compressions[b]
                if center not in compMap:
                    compMap[center] = newCluster
                    newCluster.append(b)
                else:
                    compMap[center].append(b)
            else:
                newCluster.append(b)
        if len(newCluster)>0:
            mergedClusters.append(newCluster)
            
    return mergedClusters

def countHomologies(clusters, numLetters, insertions):
    homs = 0
    for cluster in clusters:
        letters = 0
        for b in cluster:
            if b not in insertions:
                letters = letters + numLetters[b]
        homs = homs + letters * (letters - 1) / 2
    return int(homs)

def compressClustersOld(context, clusters, compressions, numLetters, threshold):
    insertions = set()
    numCells = len(context.unalignedSequences) * len(clusters)
    heap = []
    
    subIdxMap = {}
    lettersMap = {i : 0 for i in range(len(clusters))}
    heapSet = set()
    
    backCompressions = {}
    for a, b in compressions.items():
        backCompressions[b] = backCompressions.get(b, set())
        backCompressions[b].add(a)
    
    for idx, cluster in enumerate(clusters):
        for b in cluster:
            if b not in insertions:
                lettersMap[idx] = lettersMap[idx] + numLetters[b]
            subIdxMap[b] = idx
        
        dest = 0
        for b in cluster:
            if b in compressions:
                dest = max(dest, subIdxMap[compressions[b]] + 1)
            if dest >= idx:
                break
        if dest < idx:
            heapq.heappush(heap, (lettersMap[idx], idx))
            heapSet.add(idx)
    
    while len(heap) > 0 and numCells > threshold:
        ltrs, idx = heapq.heappop(heap)
        heapSet.remove(idx)
        
        destMap = {}
        move = True
        for b in clusters[idx]:
            nbr = 0
            if b in compressions:
                nbr = subIdxMap[compressions[b]] + 1
            while len(clusters[nbr]) == 0 and nbr < idx:
                nbr = nbr + 1
            destMap[b] = nbr
            if nbr >= idx:
                move = False
                break
        
        if move:
            numCells = numCells - len(context.unalignedSequences)
            for b in clusters[idx]:
                insertions.add(b)
                dest = destMap[b]
                clusters[dest].append(b)
                subIdxMap[b] = dest
                
                if b in backCompressions:
                    for nxt in backCompressions[b]:
                        nbr = subIdxMap[nxt]
                        if len(clusters[nbr]) > 0 and nbr not in heapSet:
                            heapq.heappush(heap, (lettersMap[nbr], nbr))
                            heapSet.add(nbr)
            clusters[idx] = []
    
    newClusters = [c for c in clusters if len(c) > 0]
    for c in newClusters:
        if len(c) == 1 and c[0] in insertions:
            insertions.remove(c[0])
    
    Configs.log("Compressed from {} clusters to {} clusters..".format(len(clusters), len(newClusters)))
    return newClusters, insertions

def compressClusters(context, clusters, compressions, numLetters, threshold):
    clusters = [set(c) for c in clusters]
    insertions = set()
    numCells = len(context.unalignedSequences) * len(clusters)
    heap = []
    
    subIdxMap = {}
    lettersMap = {i : 0 for i in range(len(clusters))}
    heapSet = set()
    
    backCompressions = {}
    for a, b in compressions.items():
        backCompressions[b] = backCompressions.get(b, set())
        backCompressions[b].add(a)
    
    for idx, cluster in enumerate(clusters):
        for b in cluster:
            if b not in insertions:
                lettersMap[idx] = lettersMap[idx] + numLetters[b]
            subIdxMap[b] = idx
        
        dest = 0
        for b in cluster:
            if b in compressions:
                dest = max(dest, subIdxMap[compressions[b]] + 1)
            if dest >= idx:
                break
        if dest < idx:
            heapq.heappush(heap, (lettersMap[idx], idx))
            heapSet.add(idx)
    
    while len(heap) > 0 and numCells > threshold:
        ltrs, idx = heapq.heappop(heap)
        heapSet.remove(idx)
        
        move = True
        dest = idx - 1
        while len(clusters[dest]) == 0:
            dest = dest - 1
        #if dest < 0:
        #    continue
        
        destMap = {}    
        for b in clusters[idx]:
            destMap[b] = dest
            curNode = b
            curDest = dest
            while curNode in compressions and subIdxMap[compressions[curNode]] == curDest:
                if compressions[curNode] not in insertions or curDest == 0:
                    move = False
                    break
                else:
                    curNode = compressions[curNode]
                    curDest = curDest - 1
                    while len(clusters[curDest]) == 0:
                        curDest = curDest - 1
                    destMap[curNode] = curDest
        
        if move:
            numCells = numCells - len(context.unalignedSequences)
            for b in clusters[idx]:
                insertions.add(b)
                if b in backCompressions:
                    for nxt in backCompressions[b]:
                        nbr = subIdxMap[nxt]
                        if len(clusters[nbr]) > 0 and nbr not in heapSet:
                            heapq.heappush(heap, (lettersMap[nbr], nbr))
                            heapSet.add(nbr)
                
            for b, dest in destMap.items():
                clusters[subIdxMap[b]].remove(b)
                clusters[dest].add(b)
                subIdxMap[b] = dest
    
    newClusters = [list(c) for c in clusters if len(c) > 0]
    for c in newClusters:
        if len(c) == 1 and c[0] in insertions:
            insertions.remove(c[0])
    
    Configs.log("Compressed from {} clusters to {} clusters..".format(len(clusters), len(newClusters)))
    return newClusters, insertions

def compressSubalignment(**kwargs):
    subalignmentPath = kwargs["subalignmentPath"]
    outputPath = kwargs["outputFile"]
    tempOutputPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    
    subsetAlign = sequenceutils.readFromFasta(subalignmentPath, removeDashes=False)
    subsetLen = len(next(iter(subsetAlign.values())).seq)
    
    numLetters = []    
    compressions = []
    lastIdx = {s : -1 for s in subsetAlign} 
    for i in range(subsetLen):
        notGaps = []
        for s,v in subsetAlign.items():
            if v.seq[i] != '-':
                notGaps.append(s) 
        
        numLetters.append(len(notGaps))
        dest = -1
        for s in notGaps:
            dest = max(dest, lastIdx[s])
            lastIdx[s] = i
        compressions.append(dest)
            
    with open(tempOutputPath, 'w') as textFile:
        textFile.write("{}\n".format(" ".join([str(c) for c in numLetters])))
        textFile.write("{}\n".format(" ".join([str(c) for c in compressions])))
    shutil.move(tempOutputPath, outputPath)

def writeUnconstrainedAlignment(context):
    graph = context.graph
    alignment = {}
    for taxon in context.unalignedSequences:            
        alignment[taxon] = sequenceutils.Sequence(taxon, ['-'] * len(graph.clusters))

    curIdxes = {taxon : 0 for taxon in context.unalignedSequences}
    for idx, cluster in enumerate(graph.clusters):
        for b in cluster:
            bsub, bpos = graph.matSubPosMap[b]
            taxon = context.subalignments[bsub][0]
            alignment[taxon].seq[idx] = context.unalignedSequences[taxon].seq[curIdxes[taxon]]
            curIdxes[taxon] = curIdxes[taxon] + 1

    for taxon in alignment:
        alignment[taxon].seq = "".join(alignment[taxon].seq)   
                            
    sequenceutils.writeFasta(alignment, context.outputFile)  
    Configs.log("Wrote final alignment to {}".format(context.outputFile)) 
