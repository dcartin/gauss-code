# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 12:13:26 2024

@author: cartin

This program computes the writhe of each graph diagram, for graphs obtained
from random Gauss codes by picking the node type randomly. Remember that the
writhe is only invariant under RII and RIII moves, but not RI moves, so it is
an invariant only of the graph diagram. With both graphs and their mirror
images, the average of the writhe is expected to be zero, so this is a good
test of the *graph* aspect of the random choices.

This version of the program finds the average writhe for all graphs with the
same number of *nodes*.
"""

from gauss_code import Symbol, Edge, Graph

import csv
import numpy as np

# Set number of trials for each choice of node number N, and max, min numbers
# of nodes to look at

numTrials = 5000
Nmin = 3
Nmax = 200

# Put results into lists, to write them to file later

avgWritheList = []
stdList = []

print('numTrials = {0}, Nmin = {1}, Nmax = {2}\n==='.format(numTrials, Nmin, Nmax))

for N in range(Nmin, Nmax + 1):
    
    # Keep track of writhe sum and number of accepted graphs
    
    writheList = []
    
    trials = 0
    counter = 0
    
    # Create a new graph for each trial
    
    while trials < numTrials:
        
        writhe = 0
        
        # Keep track of how many attempted graphs, to compute rejection rate
        
        counter += 1

        # Initialize with Gauss code aa
        
        G = Graph()
        G.initEdges()
    
        # Add additional nodes
    
        for iii in range(2, N + 1):
    
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
        
        # Find writhe from DT sequence
        
        [DT_seq, f_list] = G.GaussToDT()
        
        for iii in range(N):
            if f_list[DT_seq[iii][0] // 2] * DT_seq[iii][2] > 0:
                writhe += 1
            elif f_list[DT_seq[iii][0] // 2] * DT_seq[iii][2] < 0:
                writhe -= 1
                
        # # Add writhe of this graph to list
        
        writheList += [writhe]
                    
    # # Find average writhe, standard deviation
    
    avgWritheList += [np.average(writheList)]
    stdList += [np.std(writheList)]
    
    print('N = {0}'.format(N))
    
# Write results to file

with open('data-writhe-V trials = {0}.csv'.format(numTrials), 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    
    # Write headers
    
    csvwriter.writerow(['N', 'avg writhe', 'std'])
    
    for iii in range(Nmax - Nmin + 1):
        csvwriter.writerow([iii + Nmin, avgWritheList[iii], stdList[iii]])