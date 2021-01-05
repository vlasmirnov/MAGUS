'''
Created on Dec 4, 2020

@author: Vlad
'''

import os
from helpers import sequenceutils
from tasks import task

class AlignmentContext:
    
    def __init__(self, **kwargs):
        self.outputFile = None
        self.workingDir = None
        self.sequencesPath = None
        self.subsetPaths = []
        self.subalignmentPaths = []
        self.backbonePaths = []
        self.guideTreePath = None
        
        self.unalignedSequences = None
        self.subsetTaxonMap = {}
        self.taxonSubsetMap = {}
        
        self.subalignmentTasks = []
        self.graph = None
        
        for attr in kwargs:
            vars(self)[attr] = kwargs.get(attr)
        
        if not os.path.exists(self.workingDir):
            os.makedirs(self.workingDir)
    
    def initializeSequences(self):
        self.unalignedSequences = {}
        for i, subsetPath in enumerate(self.subsetPaths):
            subset = sequenceutils.readFromFasta(subsetPath, removeDashes=True)
            self.unalignedSequences.update(subset)
            self.subsetTaxonMap[i] = []
            for taxon in subset:
                self.taxonSubsetMap[taxon] = i
                self.subsetTaxonMap[i].append(taxon)
    
    def awaitSubalignments(self):
        task.awaitTasks(self.subalignmentTasks)