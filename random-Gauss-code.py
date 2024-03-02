# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 09:57:40 2024

@author: cartin

This program uses the Graph class to find a random graph with a given number of
nodes.
"""

from gauss_code import Graph

# Choose number of nodes for final graph

numNodes = 4

# Initialize with Gauss code aa

G = Graph()
G.initEdges()

# Add additional nodes

for iii in range(2, numNodes + 1):

    [edge1, edge2, face] = G.getEdgePair()
    
    if edge1 == edge2:
        G.insertSelfLoop(edge_first = edge1, face_num = face)
        
    else:
        G.insertSymbol(edge_left = edge1, edge_right = edge2, face_num = face)
        
# Print Gauss code for final graph
    
print('Gauss code:', G.lowestCode())
print('edge pair seq:', G.edge_pair_list)