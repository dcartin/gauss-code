# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 11:43:04 2024

@author: cartin

Here, we find numTrials Gauss codes with numNodes, and create random 4-regular
knotted graphs by choosing vertices and crossings with probability 1/2. Then,
for each graph obtained, we find the shortest path (length *only* in number of
edges, since there is as yet, no geometric information) between all vertices,
and find the average of all these lengths. For each graph, the number of nodes,
vertices, and the average geodesic path length are written to file.
"""

from gauss_code import Symbol, Edge, Graph

import csv

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path

# Set number of trials, and total number of nodes for each graph

numNodes = 300
numTrials = 10000

# Put results in dicts for writing to file later

avgDistDict = {}
numTrialDict = {}

# For each trial, attempt to create a graph; reject the graph if it has no
# nodes. The probability for this is 2^(-N) for N the number of nodes, so it
# is unlikely there will be a rejection beyond N ~ 5 or so.

trials = 0
counter = 0

while trials < numTrials:
    
    # Keep track of how many attempted graphs, to compute rejection rate
    
    counter += 1

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
            
    # Create types for each node; reject graphs with no vertices and start
    # from beginning of process

    G.assignNodeType()
    
    # If G has zero vertices, reject the graph and start over; otherwise,
    # record the trial, and find the stats of the graph
    
    if G.num_vert == 0:
        continue
    
    trials += 1
    
    # Find adjacency matrix, then use shortest path algorithm to find geodesic
    # distances between each vertex. There are N(N - 1)/2 total distances, and
    # each appears twice in the symmetric matrix, so divide the sum of all
    # distance matrix entries by N(N - 1) to get the average distance.
    
    M = G.adjacencyMatrix()
    numVert = len(M)
    
    distMatrix = shortest_path(csgraph = csr_matrix(M), directed = False, \
                               return_predecessors = False, unweighted = True)
    avgDist = distMatrix.sum() / (numVert * (numVert - 1))
    
    # Record results
    
    avgDistDict[numVert] = avgDistDict.get(numVert, 0) + avgDist
    numTrialDict[numVert] = numTrialDict.get(numVert, 0) + 1
    
# Write results to file

with open('data-avg-dist-N={0}.csv'.format(numNodes), 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    
    # Write headers
    
    csvwriter.writerow(['V', 'num graphs', 'avg dist'])
    
    # Write data
    
    keyList = sorted(avgDistDict.keys())
    for key in keyList:
        csvwriter.writerow([key, numTrialDict[key], avgDistDict[key] / numTrialDict[key]])