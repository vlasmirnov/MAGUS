'''
Created on Feb 4, 2020

@author: Vlad
'''

import numpy 
import os
import time
import argparse
import heapq

import external_tools
import configuration
import sequenceutils


    
def mergeAlignments(dirPath, sequencesPath, subsetPaths, guideTreePath, outputPath, mafftRuns, mafftSize):    
    unalignment = sequenceutils.readFromFasta(sequencesPath, removeDashes = True)
    #subsets = [sequenceutils.readFromFasta(filePath) for filePath in subsetPaths]
    
    subsets = []
    for filePath in subsetPaths:
        baseName = os.path.basename(filePath)
        cleanFilePath = os.path.join(dirPath, "clean_{}".format(baseName))
        sequenceutils.cleanGapColumns(filePath, cleanFilePath)
        subsets.append(sequenceutils.readFromFasta(cleanFilePath))
    
    s = len(subsets)
    subsetLengths = [len(next(iter(subsets[i].values())).seq) for i in range(s)]
    
    #matrixPath = os.path.join(dirPath, "matrix_{}".format(os.path.basename(outputPath)))
    #clusterPath = os.path.join(dirPath, "mcl_clusters_{}".format(os.path.basename(outputPath)))
    #timingsPath = os.path.join(dirPath, "timings_{}".format(os.path.basename(outputPath)))
    
    matrixPath = os.path.join(dirPath, "graph.txt")
    clusterPath = os.path.join(dirPath, "clusters.txt")
    timingsPath = os.path.join(dirPath, "timings.txt")
    
    matrixSize = sum(subsetLengths)
    matSubPosMap = getMatSubPosMap(matrixSize, s, subsetLengths)
    
    time1 = time.time() 
    if os.path.exists(matrixPath):
        thematrix = readMatrixFromFile(matrixSize, matrixPath)
    else:
        thematrix = buildMatrix(dirPath, unalignment, subsets, None, mafftRuns, mafftSize)
        writeMatrixToFile(thematrix, matrixPath)
    
    time2 = time.time()
    with open(timingsPath, 'a') as timeFile:
        timeFile.write("{} {}\n".format("Compiling backbone alignments and graph: ", time2-time1))
    
    if not os.path.exists(clusterPath):
        external_tools.runMcl(matrixPath, 4, dirPath, clusterPath)
        
    time3 = time.time()  
    with open(timingsPath, 'a') as timeFile:
        timeFile.write("{} {}\n".format("Running MCL: ", time3-time2))
    
    clusters = readClustersFromFile(clusterPath)
    orderedClusters, matSubPosMap = processClusters(matSubPosMap, clusters, thematrix)
    
    time4 = time.time()
    with open(timingsPath, 'a') as timeFile:
        timeFile.write("{} {}\n".format("Ordering clusters: ", time4-time3))
    
    finalAlignment = clustersToAlignment(orderedClusters, matSubPosMap, subsets, subsetLengths)
    #finalAlignmentFile = os.path.join(dirPath, "final_alignment_{}.txt".format(mafftRuns))
    finalAlignmentFile = os.path.join(dirPath, outputPath)
    sequenceutils.writeFasta(finalAlignment, finalAlignmentFile)

def writeMatrixToFile(thematrix, outPath):
    with open(outPath, 'w') as textFile:
        for i in range(len(thematrix)):
            for k in thematrix[i]:
                textFile.write("{} {} {}\n".format(i, k, thematrix[i][k]))

def readMatrixFromFile(matrixSize, matrixPath):
    thematrix = [{} for i in range(matrixSize)]
    with open(matrixPath) as f:
        for line in f:
            tokens = [int(token) for token in line.strip().split()]
            thematrix[tokens[0]][tokens[1]] = tokens[2]
    print("read the matrix..")
    return thematrix

def readClustersFromFile(clusterPath):
    clusters = {}
    with open(clusterPath) as f:
        for line in f:
            tokens = [int(token) for token in line.strip().split()]
            if len(tokens) > 1:
                clusters[tokens[0]] = tokens 
    print("Found {} clusters..".format(len(clusters)))
    return clusters           


def buildMatrix(dirPath, unalignment, subsets, weights, mafftRuns, mafftSize):       
    s = len(subsets)
    subsetMatrixIdx = [0] * s
    subsetLengths = [len(next(iter(subsets[i].values())).seq) for i in range(s)]
    
    for k in range(1, len(subsets)):        
        subsetMatrixIdx[k] = subsetMatrixIdx[k-1] + subsetLengths[k-1]
    
    matrixSize = sum(subsetLengths)
    thematrix = [{} for i in range(matrixSize)]
    
    minSubsetLength = min([len(subset) for subset in subsets])
    numTaxa = max(1, min(minSubsetLength, int(mafftSize/len(subsets))))
        
    for n in range(1, 1+mafftRuns):    
        joinSetFile = os.path.join(dirPath, "joinset_{}_unalign.txt".format(n))
        joinSetSubMsaTable = os.path.join(dirPath, "joinset_{}_table.txt".format(n))
        mafftResultFile = os.path.join(dirPath, "joinset_{}_mafft.txt".format(n))
        if not os.path.exists(mafftResultFile):
            joinAlign = {}
            joinTaxa = []     
            taxnum = 0
            with open(joinSetSubMsaTable, 'w') as textFile:  
                for subset in subsets:
                    taxa = list(subset.keys())
                    numpy.random.shuffle(taxa)
                    
                    for taxon in taxa[:numTaxa]:
                        taxnum = taxnum + 1
                        joinTaxa.append(taxon)
                        joinAlign[taxon] = subset[taxon]
                        textFile.write(str(taxnum) + ' ')
                    textFile.write('\n')      
               
            sequenceutils.writeFasta(joinAlign, joinSetFile, joinTaxa)
            
            external_tools.runMafft(joinSetFile, None, dirPath, mafftResultFile, 8)
        
        joinsetAlign = sequenceutils.readFromFasta(mafftResultFile)        
        alignvector, rowTaxon = joinsetToAlignVector(unalignment, joinsetAlign, subsets, subsetMatrixIdx)
        lenAlignment = len(next(iter(joinsetAlign.values())).seq) 
        
        print("Feeding backbone alignment {} to the matrix..".format(n))
        
        for l in range(lenAlignment):  
            for k1 in range(len(joinsetAlign)):
                if alignvector[k1, l] != -1:
                    for i1 in range(len(joinsetAlign)):
                        if alignvector[i1, l] != -1:
                            a, b = alignvector[k1, l], alignvector[i1, l]
                            if weights is None:
                                thematrix[a][b] = thematrix[a].get(b,0) + 1
                            else:
                                taxon1, taxon2 = rowTaxon[k1], rowTaxon[i1]
                                thematrix[a][b] = thematrix[a].get(b,0) + 1 * weights[taxon1, taxon2] 
                        
    print("compiled the matrix..")
    return thematrix
     

def joinsetToAlignVector(unalignment, joinsetAlign, subsets, subsetMatrixIdx):
    lenAlignment = len(next(iter(joinsetAlign.values())).seq)    
    alignvector = numpy.full((len(joinsetAlign), lenAlignment), -1)
    t = 0
    
    rowTaxon = [0] * len(joinsetAlign)
    for k,subset in enumerate(subsets):
        
        for taxon in subset:        
            if taxon in joinsetAlign:
                unalignedseq = unalignment[taxon].seq
                subsetseq = subset[taxon].seq
                i = 0
                posarray = [0] * len(unalignedseq)
                for n in range(len(subsetseq)):
                    if subsetseq[n] == unalignedseq[i]:
                        posarray[i] = n 
                        i = i + 1
                        if i == len(unalignedseq):
                            break
                
                joinsetseq = joinsetAlign[taxon].seq
                i = 0
                for n in range(len(joinsetseq)):
                    if i < len(unalignedseq) and joinsetseq[n] == unalignedseq[i]:
                        alignvector[t, n] = int(subsetMatrixIdx[k] + posarray[i])
                        i = i + 1
                #print(printString)
                rowTaxon[t] = taxon
                t = t + 1
            
    return alignvector, rowTaxon

def getMatSubPosMap(matrixSize, s, subsetLengths):
    matSubPosMap = [0] * matrixSize
    i = 0
    for k in range(s):
        for j in range(subsetLengths[k]):
            matSubPosMap[i] = (k, j)
            i = i + 1
    return matSubPosMap


def processClusters(matSubPosMap, clusters, thematrix):
    uniqueClusters = set()
    newclusters = {}
    for a in clusters:
        clusters[a].sort()
        clustr = tuple(clusters[a])
        if clustr not in uniqueClusters:
            uniqueClusters.add(clustr)
            newclusters[a] = clusters[a]
    clusters = newclusters
    print("Found {} unique clusters..".format(len(clusters)))
    
    redundantCols = {}
    redundantRows = {}
    elementScores = {}
    for a in clusters:
        for b in clusters[a]:
            bsub, bpos = matSubPosMap[b] 
            redundantCols[a,bsub] = redundantCols.get((a, bsub), []) + [(a, b)] 
            redundantRows[b] = redundantRows.get(b, []) + [(a,b)]
            
            thesum = 0
            for c in clusters[a]:
                if b != c:
                    csub, cpos = matSubPosMap[c]
                    if bsub != csub:
                        thesum = thesum + thematrix[b].get(c,0)
            elementScores[a, b] = thesum
    
    proburemsCols = [(a,b) for a,b in redundantCols if len(redundantCols[a,b]) > 1]
    proburemsCols.sort(key= lambda a: min([elementScores[b, c] for b, c in redundantCols[a]]))
    proburemsRows = [a for a in redundantRows if len(redundantRows[a]) > 1]
    proburemsRows.sort(key= lambda a: min([elementScores[b, c] for b, c in redundantRows[a]]))
    
    sortedScores = list(elementScores.keys())
    sortedScores.sort(key = lambda x : elementScores[x])
    
    print("Found {} and {} problems rows and cols..".format(len(proburemsRows), len(proburemsCols)))
    for a,b in sortedScores:
        bsub, bpos = matSubPosMap[b] 
        if len(redundantCols[a, bsub]) > 1:
            clusters[a].remove(b)
            #print("removed {} from cluster {} with score {}..".format(matSubPosMap[b], matSubPosMap[a], elementScores[a,b]))
            redundantCols[a, bsub].remove((a,b))
            redundantRows[b].remove((a,b))
        elif len(redundantRows[b]) > 1:
            clusters[a].remove(b)
            #print("removed {} from cluster {} with score {}..".format(matSubPosMap[b], matSubPosMap[a], elementScores[a,b]))
            redundantCols[a, bsub].remove((a,b))
            redundantRows[b].remove((a,b))
                
    
    proburemsCols = [(a,b) for a,b in redundantCols if len(redundantCols[a,b]) > 1]
    proburemsCols.sort(key= lambda a: min([elementScores[b, c] for b, c in redundantCols[a]]))
    proburemsRows = [a for a in redundantRows if len(redundantRows[a]) > 1]
    proburemsRows.sort(key= lambda a: min([elementScores[b, c] for b, c in redundantRows[a]]))
    print("But now {} and {} problems rows and cols..".format(len(proburemsRows), len(proburemsCols)))
    
    clusters = [clusters[a] for a in clusters if len(clusters[a]) > 1]
    print("Found {} purified clusters..".format(len(clusters)))

    orderedClusters = orderClusters(clusters, matSubPosMap)
            
    print("")
    print("Compiled {} sorted clusters..".format(len(orderedClusters)))
    #for n, cluster in enumerate(orderedClusters):
    #    print(n, [matSubPosMap[b] for b in cluster])

    return orderedClusters, matSubPosMap

def orderClusters(clusters, matSubPosMap):
    subsetClusters = {}
    clusterPositions = {}
    queueIdxs = {}
    clusterBreaks = {}
    maxFrontier = {}
    visitedBreaks = {}
    visitedStates = set()
    maximalCut = {}
    trashCounter = 0
    heapLimit = 5000
    aggression = 1.0
    lastFrontierState = None
    branchingLimit = float('inf')
    totalPairs = 0

    
    for a,cluster in enumerate(clusters):
        for b in cluster:
            bsub, bpos = matSubPosMap[b] 
            subsetClusters[bsub] = subsetClusters.get(bsub, []) + [(a, bpos)]  
            clusterPositions[a] = {}
              
            visitedBreaks[a] = set()
            totalPairs = totalPairs + len(cluster)*(len(cluster)-1)/2
                
    
    for asub in subsetClusters:
        subsetClusters[asub].sort(key = lambda c: c[1])
        queueIdxs[asub] = 0

        maxFrontier[asub] = -1
        maximalCut[asub] = -1
        
        for i in range(len(subsetClusters[asub])):            
            a = subsetClusters[asub][i][0]
            clusterPositions[a][asub] = i
    
    heap = []
    
    startState = (0, 0, len(clusters), totalPairs, trashCounter, queueIdxs, clusterBreaks, maximalCut, [], True)
    #startState = (0, 0, totalPairs, trashCounter, queueIdxs, clusterBreaks, maximalCut, [], True)
    startState = developState(startState, aggression, branchingLimit, 0, clusters, subsetClusters, clusterPositions, matSubPosMap)
    heapq.heappush(heap, startState)
    
    while len(heap) > 0:
        cleared = False
        if len(heap) > heapLimit:
            print("Heap limit {} reached.. Truncating heap to last frontier".format(heapLimit)) 
            if aggression == 1:
                aggression = 1.2
                print("Increasing aggression to {}..".format(aggression))
            elif aggression < 8:
                aggression = int(aggression) * 2
                print("Increasing aggression to {}..".format(aggression))
            else:
                print("Setting strategy to fully greedy..")
                branchingLimit = 1
                aggression = 1    
                
            heap.clear()
            visitedStates = set()
            lastFrontierState = developState(lastFrontierState, aggression, branchingLimit, 0, clusters, subsetClusters, clusterPositions, matSubPosMap)
            heapq.heappush(heap, lastFrontierState)
            cleared = True

        state = heapq.heappop(heap)
        heuristic, numOrdered, numLeft, pairsLeft, trash, queueIdxs, clusterBreaks, maximalCut, newClusterBreaks, safeFrontier = state
   
        if len(newClusterBreaks) == 0:
            break
        else:
            stateKey = tuple([queueIdxs[sub] for sub in subsetClusters])
            if stateKey in visitedStates:
                #print("alrdy visited..")
                continue
            else:
                visitedStates.add(stateKey)
            

            newSubFrontier = True
            for asub in queueIdxs:
                if queueIdxs[asub] <= maxFrontier[asub]:
                    newSubFrontier = False
                    break
                
            if newSubFrontier:
                maxFrontier = queueIdxs
                print("New frontier")
                print(maxFrontier)
                lastFrontierState = state
                
                branchingLimit = float('inf')
                #aggression = 1.0
                #heap.clear()
                #visitedStates = set()
            
            if safeFrontier and not cleared:
                print("Safe frontier.. dumping {} from heap and resetting aggression..".format(len(heap)))
                lastFrontierState = state
                
                heap.clear()
                visitedStates = set()
                aggression = 1.0
                branchingLimit = float('inf')
                #heapq.heappush(heap, lastFrontierState)
                #continue
                
            
            nextStates = []
            for a, goodSide, badSide, crossedClusters in newClusterBreaks:
                
                g,b = len(goodSide), len(badSide)
                pairsDiff = g*(g-1)/2 + b*(b-1)/2 - (g+b)*(g+b-1)/2
                
                trashCounter = trashCounter + 1
                queueIdxsCopy = dict(queueIdxs)
                clusterBreaksCopy = dict(clusterBreaks)
                maximalCutCopy = dict(maximalCut)
                for b in goodSide:
                    bsub, bpos = matSubPosMap[b]
                    clusterBreaksCopy[a, bsub] = goodSide
                    maximalCutCopy[bsub] = max(maximalCutCopy[bsub], clusterPositions[a][bsub])
                    
                for b in badSide:
                    bsub, bpos = matSubPosMap[b]
                    clusterBreaksCopy[a, bsub] = badSide
                    maximalCutCopy[bsub] = max(maximalCutCopy[bsub], clusterPositions[a][bsub])
                
                nextState = (0, numOrdered, numLeft + 1, pairsLeft + pairsDiff, trashCounter, queueIdxsCopy, clusterBreaksCopy, maximalCutCopy, [], False)
                nextState = developState(nextState, aggression, branchingLimit, len(crossedClusters), clusters, subsetClusters, clusterPositions, matSubPosMap)

                nextStates.append(nextState)
                
            if branchingLimit < len(nextStates):
                nextStates.sort(key=lambda x : x[0], reverse=False)
                
            for num, nextState in enumerate(nextStates):
                if num >= branchingLimit:
                    break
                heapq.heappush(heap, nextState)
                
    queueIdxs = {}
    for asub in subsetClusters:
        queueIdxs[asub] = 0
        
    orderedClusters = []
    foundGood = True
    while foundGood:
        foundGood = False
        for asub in queueIdxs:
            good = True
            idx = queueIdxs[asub]
            if idx == len(subsetClusters[asub]):
                continue
            a, pos = subsetClusters[asub][idx]
            
            if (a,asub) in clusterBreaks:
                cluster = clusterBreaks[a,asub]
            else:
                cluster = clusters[a]
            
            for b in cluster:
                bsub, bpos = matSubPosMap[b]
                if clusterPositions[a][bsub] != queueIdxs[bsub]:
                    good = False
                    break
                
            if good:
                orderedClusters.append(cluster)
                for b in cluster:
                    bsub, bpos = matSubPosMap[b]
                    queueIdxs[bsub] = clusterPositions[a][bsub] + 1
                foundGood = True
                break
                    
    return orderedClusters

def developState(state, aggression, branchingLimit, crossed, clusters, subsetClusters, clusterPositions, matSubPosMap):
    heuristic, numOrdered, numLeft, pairsLeft, trash, queueIdxs, clusterBreaks, maximalCut, newClusterBreaks, safeFrontier  = state
    
    foundGood = True
    while foundGood:
        foundGood = False
        newClusterBreaks = []
        visited = set()
        safeFrontier = True
        
        for asub in queueIdxs:
            idx = queueIdxs[asub]
            if idx <= maximalCut[asub]:
                safeFrontier = False

            if idx == len(subsetClusters[asub]):
                continue
            a, pos = subsetClusters[asub][idx]
            if (a,asub) in visited:
                continue
            
            if (a,asub) in clusterBreaks:
                cluster = clusterBreaks[a,asub]
            else:
                cluster = clusters[a]
            
            goodSide, badSide, crossedClusters = [], [], set()
            for b in cluster:
                bsub, bpos = matSubPosMap[b]
                visited.add((a, bsub))
                bidx = clusterPositions[a][bsub]
                diff = bidx - queueIdxs[bsub]
                if diff == 0:
                    goodSide.append(b)
                else:
                    badSide.append(b)
                                   
            if len(badSide) == 0:
                for b in cluster:
                    bsub, bpos = matSubPosMap[b]
                    queueIdxs[bsub] = clusterPositions[a][bsub] + 1
                numOrdered = numOrdered + 1
                numLeft = numLeft - 1
                foundGood = True
                break
            else:
                newClusterBreaks.append((a, goodSide, badSide, crossedClusters))

    
    if branchingLimit == 1:
        for a, goodSide, badSide, crossedClusters in newClusterBreaks:
            goodSub = set()
            for b in goodSide:
                bsub, bpos = matSubPosMap[b] 
                goodSub.add(bsub)
                
            for b in badSide:    
                bsub, bpos = matSubPosMap[b]   
                for i in range(queueIdxs[bsub], clusterPositions[a][bsub]):
                    c, posc = subsetClusters[bsub][i]
                    if (c,bsub) in clusterBreaks:
                        otherCluster = clusterBreaks[c,bsub]
                    else:
                        otherCluster = clusters[c]
                    
                    for csite in otherCluster:
                        csub, cpos = matSubPosMap[csite]
                        if csub in goodSub and clusterPositions[c][csub] > clusterPositions[a][csub]:
                            crossedClusters.add(c)
                            break  
                          
                
    if safeFrontier or len(newClusterBreaks) == 0:
        heuristic = (numLeft + numOrdered, -numOrdered, -crossed, -pairsLeft)
    else:
        heuristic = (aggression * numLeft + numOrdered, -numOrdered, -crossed, -pairsLeft)
    state = (heuristic, numOrdered, numLeft, pairsLeft, trash, queueIdxs, clusterBreaks, maximalCut, newClusterBreaks, safeFrontier)
    return state

def clustersToAlignment(orderedClusters, matSubPosMap, subsets, subsetLengths):
    finalAlignment = {}
    for subset in subsets:
        for taxon in subset:
            finalAlignment[taxon] = sequenceutils.Sequence(taxon, "")
    
    lastIdx = {}
    for cluster in orderedClusters:
        for b in cluster:
            bsub, bpos = matSubPosMap[b]
            bposlast = lastIdx.get(bsub, -1)
            
            for p in range(bposlast + 1, bpos):
                for i in range(len(subsets)):
                    for taxon in subsets[i]:
                        if i == bsub:
                            finalAlignment[taxon].seq = finalAlignment[taxon].seq + subsets[bsub][taxon].seq[p]
                        else:
                            finalAlignment[taxon].seq = finalAlignment[taxon].seq + '-'
            
            lastIdx[bsub] = bpos
        
        clusterSets = set()    
        for b in cluster:
            bsub, bpos = matSubPosMap[b]
            clusterSets.add(bsub)
            for taxon in subsets[bsub]:
                finalAlignment[taxon].seq = finalAlignment[taxon].seq + subsets[bsub][taxon].seq[bpos]
            #lastIdx[bsub] = lastIdx[bsub] + 1
        
        for i in range(len(subsets)):
            if i not in clusterSets:
                for taxon in subsets[i]:
                    finalAlignment[taxon].seq = finalAlignment[taxon].seq + '-'
            
    
    for bsub in range(len(subsets)):
        bposlast = lastIdx.get(bsub, -1)
                    
        for p in range(bposlast + 1, subsetLengths[bsub]):
            for i in range(len(subsets)):
                for taxon in subsets[i]:
                    if i == bsub:
                        finalAlignment[taxon].seq = finalAlignment[taxon].seq + subsets[bsub][taxon].seq[p]
                    else:
                        finalAlignment[taxon].seq = finalAlignment[taxon].seq + '-'        
    return finalAlignment

def main(args):   
    startTime = time.time()
    
    configuration.buildConfigs(args)
    workingDir = args.directory   
    sequencesPath = args.sequences
    subsetPaths = args.subalignments
    guideTree = args.guidetree
    outputPath = args.output
    mafftRuns = args.mafftruns
    mafftSize = args.mafftsize
    
    mergeAlignments(workingDir, sequencesPath, subsetPaths, guideTree, outputPath, mafftRuns, mafftSize)
    
    endTime = time.time()
    print("Finished aligning in {} seconds..".format(endTime-startTime))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str,
                        help="Path to working directory", required=True)
    
    parser.add_argument("-s", "--sequences", type=str,
                        help="Path to unaligned sequences", required=True)

    parser.add_argument("-a", "--subalignments", type=str, nargs="+",
                        help="Paths to subset alignment files", required=True)

    parser.add_argument("-o", "--output", type=str,
                        help="Output alignment path",
                        required=True)
    
    parser.add_argument("-t", "--guidetree", type=str,
                        help="guide tree for merge weights",
                        required=False, default=None)
    
    parser.add_argument("-r", "--mafftruns", type=int,
                        help="Number of MAFFT runs", required=False, default=10)
    
    parser.add_argument("-m", "--mafftsize", type=int,
                        help="Maximum size of MAFFT alignments", required=False, default=200)

    main(parser.parse_args())