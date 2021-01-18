'''
Created on Jan 6, 2021

@author: Vlad
'''

import os

from configuration import Configs
from helpers import sequenceutils
from tasks import task
import shutil

def writeAlignment(context):
    context.graph.addSingletonClusters()
    if Configs.constrain:
        writePackedAlignment(context)    
        writeUnpackedAlignment(context)
    else:
        writeUnconstrainedAlignment(context)

def writePackedAlignment(context):
    graph = context.graph
    packedAlignment = {}
    for path in context.subalignmentPaths:            
        packedAlignment[path] = sequenceutils.Sequence(path, ['-'] * len(graph.clusters))

    for idx, cluster in enumerate(graph.clusters):
        for b in cluster:
            bsub, bpos = graph.matSubPosMap[b]
            packedAlignment[context.subalignmentPaths[bsub]].seq[idx] = 'X'

    for s in packedAlignment:
        packedAlignment[s].seq = "".join(packedAlignment[s].seq)   
                            
    sequenceutils.writeFasta(packedAlignment, graph.packedAlignmentPath)  
    Configs.log("Wrote uncompressed packed subset alignment to {}".format(graph.packedAlignmentPath))
    if not os.path.exists(graph.alignmentMaskPath):  
        buildAlignmentMask(context)
    
def buildAlignmentMask(context):
    graph = context.graph
    numCells = len(context.unalignedSequences) * len(graph.clusters)
    Configs.log("Final alignment will have {} cells..".format(numCells))
    Configs.log("Building alignment mask..")
    if numCells > 5e10:
        Configs.log("Removing columns with more than 99% gaps..")
        numSeq = len(context.unalignedSequences)
        portion = 0.99
        mask = [0] * len(graph.clusters) 
        gapCounts = buildGapCounts(context)
        
        for n, cluster in enumerate(graph.clusters):
            nonGaps = 0
            for b in cluster:
                bsub, bpos = graph.matSubPosMap[b]
                nonGaps = nonGaps + len(context.subalignments[bsub]) - gapCounts[context.subalignmentPaths[bsub]][bpos]
                if nonGaps >= (1-portion) * numSeq:
                    mask[n] = 1
                    break       
    else:
        mask = [1] * len(graph.clusters)
        
    maskString = "".join([str(c) for c in mask])
    with open(graph.alignmentMaskPath, 'w') as u:
        u.write("{}\n".format(maskString))
    
def buildGapCounts(context):
    countingTasks = []
    for subalignPath in context.subalignmentPaths:
        countsPath = os.path.join(context.graph.workingDir, "gap_counts_{}".format(os.path.basename(subalignPath)))
        args = {"alignFile" : subalignPath, "outputFile" : countsPath}
        countingTask = task.Task(taskType = "recordGapCounts", outputFile = args["outputFile"], taskArgs = args)
        countingTasks.append(countingTask)
    
    if len(context.unalignedSequences) >= 10000:
        task.submitTasks(countingTasks)
    else:
        for countingTask in countingTasks:
            countingTask.run()
            countingTask.isFinished = True
    
    gapCounts = {}
    for countingTask in task.asCompleted(countingTasks):
        with open(countingTask.outputFile) as f:
            countLine = [int(c) for c in f.readline().strip().split()]
            gapCounts[countingTask.taskArgs["alignFile"]] = countLine
        os.remove(countingTask.outputFile)    
    return gapCounts

def writeUnpackedAlignment(context):
    graph = context.graph
    filePath = context.outputFile
    tempFile = os.path.join(os.path.dirname(filePath), "temp_{}".format(os.path.basename(filePath)))
    if os.path.exists(tempFile):
        os.remove(tempFile)
        
    tempMaskFile = os.path.join(graph.workingDir, "temp_{}".format(os.path.basename(graph.maskedSequencesPath)))
    if os.path.exists(tempMaskFile):
        os.remove(tempMaskFile)
        
    Configs.log("Assembling final alignment in {}".format(filePath))
    inducedSubalignTasks = []
    for subalignPath in context.subalignmentPaths:
        inducedAlignPath = os.path.join(graph.workingDir, "induced_{}".format(os.path.basename(subalignPath)))
        maskedSeqPath = os.path.join(graph.workingDir, "masked_{}".format(os.path.basename(subalignPath)))
        args = {"packedFilePath" : graph.packedAlignmentPath, "subalignmentPath" : subalignPath, "outputFile" : inducedAlignPath, 
                "alignmentMaskPath" : graph.alignmentMaskPath, "maskedSequencesPath" : maskedSeqPath}
        inducedTask = task.Task(taskType = "buildInducedSubalignment", outputFile = args["outputFile"], taskArgs = args)
        inducedSubalignTasks.append(inducedTask)
    
    if len(context.unalignedSequences) >= 10000:
        task.submitTasks(inducedSubalignTasks)
    else:
        for inducedTask in inducedSubalignTasks:
            inducedTask.run()
            inducedTask.isFinished = True
            
    for inducedTask in task.asCompleted(inducedSubalignTasks):
        inducedAlign = sequenceutils.readFromFasta(inducedTask.outputFile, removeDashes=False)
        Configs.log("Appending induced alignment, {} sequences of length {}..".format(len(inducedAlign), len(next(iter(inducedAlign.values())).seq)))
        sequenceutils.writeFasta(inducedAlign, tempFile, append = True)   
        
        maskedSeq = sequenceutils.readFromFasta(inducedTask.taskArgs["maskedSequencesPath"])
        sequenceutils.writeFasta(maskedSeq, tempMaskFile, append = True)
          
        os.remove(inducedTask.outputFile)
        os.remove(inducedTask.taskArgs["maskedSequencesPath"])
    shutil.move(tempMaskFile, graph.maskedSequencesPath)
    shutil.move(tempFile, filePath)
    Configs.log("Wrote final alignment to {}".format(filePath))    

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
    
def recordGapCounts(**kwargs):
    counts = sequenceutils.countGaps(kwargs["alignFile"])
    countString = " ".join([str(c) for c in counts])
    with open(kwargs["outputFile"], 'w') as u:
        u.write("{}\n".format(countString))
    
def buildInducedSubalignment(**kwargs):
    packedFilePath = kwargs["packedFilePath"]
    subalignmentPath = kwargs["subalignmentPath"]
    alignmentMaskPath = kwargs["alignmentMaskPath"]
    maskedSequencesPath = kwargs["maskedSequencesPath"]
    inducedAlignPath = kwargs["outputFile"]
    
    tempMaskedSequencesPath = os.path.join(os.path.dirname(maskedSequencesPath), "temp_{}".format(os.path.basename(maskedSequencesPath)))
    tempInducedAlignPath = os.path.join(os.path.dirname(inducedAlignPath), "temp_{}".format(os.path.basename(inducedAlignPath)))
    
    packedAlign = sequenceutils.readFromFasta(packedFilePath, removeDashes=False)
    subsetAlign = sequenceutils.readFromFasta(subalignmentPath, removeDashes=False)
    mask = packedAlign[subalignmentPath].seq
    inducedAlign = {taxon : [] for taxon in subsetAlign}
    maskedAlignment = {taxon : [] for taxon in subsetAlign}
    with open(alignmentMaskPath) as f:
        inclusionMask = [int(c) for c in f.readline().strip()]
        
    j = 0
    for i in range(len(mask)):
        if mask[i] == '-' and inclusionMask[i] == 1:
            for taxon in subsetAlign:
                inducedAlign[taxon].append('-')
        elif mask[i] != '-':
            for taxon in subsetAlign:
                letter = subsetAlign[taxon].seq[j]
                if inclusionMask[i] == 1:
                    inducedAlign[taxon].append(letter)
                if letter != '-':
                    maskedAlignment[taxon].append(letter.upper() if inclusionMask[i] == 1 else letter.lower())
            j = j + 1
    
    for s in maskedAlignment:
        maskedAlignment[s] = sequenceutils.Sequence(s, "".join(maskedAlignment[s]))   
    sequenceutils.writeFasta(maskedAlignment, tempMaskedSequencesPath)
    shutil.move(tempMaskedSequencesPath, maskedSequencesPath)
            
    for s in inducedAlign:
        inducedAlign[s] = sequenceutils.Sequence(s, "".join(inducedAlign[s]))
    sequenceutils.writeFasta(inducedAlign, tempInducedAlignPath)
    shutil.move(tempInducedAlignPath, inducedAlignPath)