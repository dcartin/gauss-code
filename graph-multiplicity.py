# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 08:33:32 2024

@author: cartin

This program creates all Gauss codes with a given number of nodes, and finds
the multiplicity of each isomorphic set of graphs.
"""

from copy import deepcopy
from gauss_code import Graph
from Graphpoly import modOrbit

# Choose number of nodes for final graph

numNodes = 5

# Initialize with Gauss code aa

G = Graph()
G.initEdges()

queue = [G]
finalList = []

# Add additional nodes

while len(queue) > 0:
    
    # Pull graph off of queue, and find all possible Gauss codes with one
    # larger number of nodes
    
    currentGraph = queue.pop()
    
    numPairs = currentGraph.numPairs()

    for pairNum in range(numPairs):
        nextGraph = deepcopy(currentGraph)
        
        [edge1, edge2, face] = nextGraph.getEdgePair(pairNum)
        
        if edge1 == edge2:
            nextGraph.insertSelfLoop(edge_first = edge1, face_num = face)
            
        else:
            nextGraph.insertSymbol(edge_left = edge1, edge_right = edge2, face_num = face)
        
        # Depending on number of nodes, add to queue or final graph
        
        if nextGraph.num_nodes == numNodes:
            finalList += [nextGraph]
        else:
            queue += [nextGraph]
            
# Find multiplicity of each resulting graph, using the DT sequence for each graph

multDict = {}

for graph in finalList:
    DTseq = modOrbit(graph.GaussToDT(f_list = False)[0], crossing = False)
    DTstring = ' '.join([repr(node[1]) for node in DTseq])

    multDict[DTstring] = multDict.get(DTstring, 0) + 1
    
# Sort as numbers, not as strings (which was necessary to use for dict)

seqList = []

for key in multDict.keys():
    seqList += [(key, [int(label) for label in key.split(' ')])]
    
seqList.sort(key = lambda seq : seq[1])

for ting in seqList:
    print(ting[0], multDict[ting[0]])
    
print('num of nodes:', numNodes)
print('num of non-iso seqs:', len(multDict))
print('num of graphs:', len(finalList))