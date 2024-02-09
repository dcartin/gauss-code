# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 08:33:32 2024

@author: cartin

This program creates all Gauss codes with a given number of nodes, and finds
the multiplicity of each isomorphic set of graphs.

09 Feb 2024: Added part to write to file the multiplicity and edge pair list
for each non-isomorphic Gauss code.
"""

from copy import deepcopy
from gauss_code import Graph
from Graphpoly import modOrbit

from time import time
start = time()

# Choose number of nodes for final graph

numNodes = 5

# Initialize with Gauss code aa

G = Graph()
G.initEdges()

currentDict = {G : 1}

# Add additional nodes

for iii in range(2, numNodes + 1):
    
    print('N =', iii)
    
    graphDict = {}
    numDict = {}
    
    # Go through all current non-isomorphic graphs, and find all graphs with
    # one greater node

    for currentGraph in currentDict.keys():
        
        numPairs = currentGraph.numPairs()
    
        for pairNum in range(numPairs):
            nextGraph = deepcopy(currentGraph)
            
            [edge1, edge2, face] = nextGraph.getEdgePair(pairNum)
            
            if edge1 == edge2:
                nextGraph.insertSelfLoop(edge_first = edge1, face_num = face)
                
            else:
                nextGraph.insertSymbol(edge_left = edge1, edge_right = edge2, face_num = face)
            
            # Add next graph to dicts only if it has a distinct lowest order
            # Gauss code; be sure to include multiplicity of previous graph
            
            lowest = nextGraph.lowestCode()
            
            if lowest not in graphDict.keys():
                graphDict[lowest] = nextGraph
                numDict[lowest] = currentDict[currentGraph]
            else:
                numDict[lowest] += currentDict[currentGraph]
            
    # Move all graphs to currentList, with correct multiplicity
            
    currentDict = {}
    
    for code in graphDict.keys():
        currentDict[graphDict[code]] = numDict[code]
        
print('---')
        
# For desired number of nodes, print out all Gauss codes, multiplicities, and
# edge pair number lists

edgePairList = [(numDict[code], graphDict[code].edge_pair_list) for code in graphDict.keys()]
edgePairList.sort(key = lambda pair : pair[1])

with open('edgePairList N = {}.txt'.format(numNodes), 'w') as f:
    
    for item in edgePairList:
        f.write(repr(item) + '\n')
        
f.close()

print('num of distinct Gauss codes:', len(graphDict.keys()), '\n---')
        
# Find multiplicity of each resulting graph, using the DT sequence for each graph

multDict = {}
pairListDict = {}

for graph in currentDict.keys():
    DTseq = modOrbit(graph.GaussToDT(f_list = False)[0], crossing = False)
    DTstring = ' '.join([repr(node[1]) for node in DTseq])

    multDict[DTstring] = multDict.get(DTstring, 0) + currentDict[graph]
    pairListDict[DTstring] = graph.edge_pair_list
    
# Sort as numbers, not as strings (which was necessary to use for dict)

seqList = []

for key in multDict.keys():
    seqList += [(key, [int(label) for label in key.split(' ')], pairListDict[key])]
    
seqList.sort(key = lambda seq : seq[1])

for seq in seqList:
    print('{0} \t {1:5d} \t {2}'.format(seq[0], multDict[seq[0]], seq[2]))
    
print('num of nodes:', numNodes)
print('num of non-iso seqs:', len(multDict))
print('num of graphs:', sum(currentDict.values()))
    
print('time:', time() - start)