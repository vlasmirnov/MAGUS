'''
Created on Apr 13, 2020

@author: Vlad
'''

import os

from configuration import Configs
from tools import external_tools


def runMlrMclClustering(graph):  
    Configs.log("Running MLR-MCL alignment graph clustering..")
    clusterPath = os.path.join(graph.workingDir, "clusters_mlr_mcl.txt")
    
    if not os.path.exists(clusterPath):
        graphPath = os.path.join(graph.workingDir, "library_graph_mlr_mcl.txt")
        if not os.path.exists(graphPath):
            writeGraphToFile(graph, graphPath)
        else:
            Configs.log("Found existing MLR-MCL graph file {}".format(graphPath))
        
        external_tools.runMlrMcl(graphPath, 30000, 0.5, 4, graph.workingDir, clusterPath)
        #external_tools.runMlrMcl(graphPath, None, None, None, graph.workingDir, clusterPath)
        #external_tools.runMlrMcl(graphPath, 1000000, None, None, graph.workingDir, clusterPath)
        #external_tools.runMlrMcl(graphPath, None, None, 8, graph.workingDir, clusterPath)
        #external_tools.runMlrMcl(graphPath, None, 0.25, None, graph.workingDir, clusterPath)
    else:
        Configs.log("Found existing MLR-MCL cluster file {}".format(clusterPath))
    graph.clusters = readClustersFromFile(clusterPath)
    
    
def writeGraphToFile(graph, filePath):
    Configs.log("Writing MLR-MCL graph file to {}".format(filePath))
    vertices, edges = 0, 0
    lines = []
    for i in range(len(graph.matrix)):
        pairs = graph.matrix[i].items()
        vertices = vertices + 1        
        edges = edges + len(pairs)
        lines.append(" ".join(["{} {}".format(a+1, b) for a, b in pairs]))      
            
    with open(filePath, 'w') as textFile:
        textFile.write("{} {} 1\n".format(vertices, int(edges/2)))
        for line in lines:
            textFile.write(line + "\n")
            
    Configs.log("Wrote graph with {} vertices and {} edges to {}".format(vertices, int(edges/2), filePath))

def readClustersFromFile(filePath):
    assignments = {}
    with open(filePath) as f:
        num = 0
        for line in f:
            cluster = int(line.strip())
            if cluster not in assignments:
                assignments[cluster] = [num]
            else:
                assignments[cluster].append(num) 
            num = num + 1
    clusters = [assignments[c] for c in range(len(assignments))]
    Configs.log("Found {} clusters..".format(len(clusters)))
    return clusters