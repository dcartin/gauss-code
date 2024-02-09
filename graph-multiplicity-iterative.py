# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 08:33:32 2024

@author: cartin

This program creates all Gauss codes with a given number of nodes, and finds
the multiplicity of each isomorphic set of graphs.

09 Feb 2024: This version uses the list of all graphs with one less node, to
build up the next node number level. The values seem to be correct when I use
all non-isomorphic orientation-signed Gauss codes, rather than their DT codes
(without orientation information).
"""

from copy import deepcopy
from gauss_code import Graph
from Graphpoly import modOrbit

from time import time
start = time()

# Open file with multiplicities and edge pair lists for one lower number of nodes

numNodes = 4            # Number of nodes to *read in*, not create
edgePairList = []

with open('edgePairList N = {}.txt'.format(numNodes), 'r') as f:
    for line in f:
        edgePairList += [eval(line.rstrip())]
        
f.close()

numNodes += 1           # Number of nodes to *create*

# Go through list, and create graphs with one greater node

graphDict = {}
numDict = {}

for (seedNum, (mult, seedList)) in enumerate(edgePairList):
    print('pair num:', seedNum)

    # Initialize with Gauss code aa
    
    G = Graph()
    G.initEdges()
    
    # Go through seedList, and create graph
    
    for edgePair in seedList:
        [edge1, edge2, face] = G.getEdgePair(edgePair)
        
        if edge1 == edge2:
            G.insertSelfLoop(edge_first = edge1, face_num = face)
            
        else:
            G.insertSymbol(edge_left = edge1, edge_right = edge2, face_num = face)
            
    # With graph created, go through all possible ways to add a node
    
    numPairs = G.numPairs()
    
    for pairNum in range(numPairs):
        nextGraph = deepcopy(G)
        
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
            numDict[lowest] = mult
        else:
            numDict[lowest] += mult
        
    # Move all graphs to currentList, with correct multiplicity
            
    currentDict = {}
    
    for code in graphDict.keys():
        currentDict[graphDict[code]] = numDict[code]
        
# For desired number of nodes, print to file all multiplicities and edge pair
# number lists, for each non-isomorphic orientation-signed Gauss code

edgePairList = [(numDict[code], graphDict[code].edge_pair_list) for code in graphDict.keys()]
edgePairList.sort(key = lambda pair : pair[1])

with open('edgePairList N = {}.txt'.format(numNodes), 'w') as f:
    for item in edgePairList:
        f.write(str(item) + '\n')
    
f.close()
        
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
    print('{0} \t {1:6d} \t {2}'.format(seq[0], multDict[seq[0]], seq[2]))
    
print('num of nodes:', numNodes)
print('num of non-iso seqs:', len(multDict))
print('num of graphs:', sum(currentDict.values()))
    
print('time:', time() - start)