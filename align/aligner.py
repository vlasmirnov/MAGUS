'''
Created on May 29, 2020

@author: Vlad
'''

import os
import time

from align.decompose import decomposer 
from align.merge.merger import mergeSubalignments
from tools import external_tools
from configuration import Configs

def alignSequences(workingDir, sequencesPath, outputPath):
    time1 = time.time()
    
    subsetPaths = decomposer.decomposeSequences(workingDir, sequencesPath)
    
    time2 = time.time()  
    Configs.log("Decomposed {} into {} subsets in {} sec..".format(sequencesPath, len(subsetPaths), time2-time1))
    Configs.log("Building {} subalignments for {}".format(len(subsetPaths), sequencesPath))
    
    subsetAlignMap = {path : os.path.join(os.path.dirname(path), "mafft_{}".format(os.path.basename(path))) for path in subsetPaths}
    external_tools.buildMafftAlignments(subsetAlignMap)
    subalignmentPaths = list(subsetAlignMap.values())
    
    time3 = time.time()  
    Configs.log("Built {} subalignments for {} in {} sec..".format(len(subalignmentPaths), sequencesPath, time3-time2))
    
    mergeSubalignments(workingDir, subalignmentPaths, outputPath)
    
    time4 = time.time()  
    Configs.log("Merged the {} subalignments for {} in {} sec..".format(len(subalignmentPaths), sequencesPath, time4-time3))