# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 10:03:35 2024

@author: cartin

This program gives the necessary classes to use node additions on Gauss codes.
It creates a class for each edge, and the graph as a whole. This is a modified
version of tbe original, which had classes for symbols and faces as well. Here,
all information is placed either in the edge class (which face the edge is part
of, and its rotation around that edge) or a dictionary of which edges correspond
to a particular face.

Notation:
    
    (1) Nodes are given by numbers +/- 1, 2, 3... The sign represents the
    orientation of the other edge crossing at the same node. If the other edge
    crosses right to left, the label is positive; otherwise it is negative.
    
    NOTE: this is the opposite of the lecture notes from Jeff Erickson, but
    seems to be a more standard nomenclature. In particular, it means the sign
    of the node matches the sign of f(i) from Dowker and Thistletwaite.
    
    (2) Edges are given by numbers 0, 1, 2, ... with no sign. The orientation
    of the edge comes from the first and last attributes, which give the
    appropriate node labels.
    
    (3) Faces are represented by unsigned numbers 1, 2, 3, ... when appearing
    in the face_left and face_right attributes of the Edge class. However, in
    faceDict and edgeNumDict for the Graph class, the numbers are signed, to
    show whether the face appears to the left of the edges listed (positive
    face label) or to the right (negative face label).
    
"""

from itertools import combinations_with_replacement, permutations
from random import randrange

#-----------------------------------------------------------------------------#

class Symbol:
    """
        A node is given by a label 1, 2, 3, ... and a sign, depending on how
        the Eulerian circuit passes through the node. When passing through the
        node, such that the other part of the circuit is right-to-left, then
        the symbol is positive; otherwise, it is negative.
    """
    
    def __init__(self, label = None, chord = None):
        
        self.label = label
        self.chord = chord
        
    def __repr__(self):
        
        if self.label > 0:
            return repr(self.label) + '+'
        
        return repr(-self.label) + '-'
        
    def linkSymbols(self, other):
        """
            Connect two symbols together to represent the same node.
        """
        
        self.chord = other
        other.chord = self
        
    def revSign(self):
        """
            When the orientation of one of the two edges passing through a
            node is reversed, this flips the sign of the two symbols for both
            edges.
        """
        
        self.label = -self.label
        self.chord.label = -self.chord.label

#-----------------------------------------------------------------------------#
        
class Edge:
    """
        An edge in the graph is represented by two symbols. The symbols are
        given in the order they are encountered along the Eulerian circuit
        defined by the Gauss code.
    """
    
    def __init__(self, first = None, last = None, face_left = None, \
                 face_right = None, prev = None, next = None, edge_num = None):
        
        # The symbols representing the nodes the edge is incident to; the
        # orientation is a -> b, where the ordered pair is given as (a, b)
        
        self.first = first
        self.last = last
        
        # The faces which the edge is part of; "left" and "right" are relative
        # to moving along the edge with its orientation. The edge will travel
        # CW around the left face, and CCW around the right.
        
        self.face_left = face_left
        self.face_right = face_right
        
        # Information about how the edge is placed in the linked list of edges
        
        self.prev = prev
        self.next = next
        
        # Keep track of edges using unique number
        
        self.edge_num = edge_num
        
    def __repr__(self):
        """
            Return edge as ordered pair, with the order giving the orientation
            of the edge first -> last.
        """
        
        # return repr(self.edge_num) + ': (' + repr(self.first) + ', ' + repr(self.last) + ')'
        return repr(self.edge_num)
    
    def __lt__(self, other):
        """
            Enable ordering of the edges, based on their label.
        """
        
        return self.edge_num < other.edge_num

#-----------------------------------------------------------------------------#
    
class Graph:
    
    def __init__(self):
        
        self.head = None
        
        # These numbers should always be related by the Euler characteristic,
        # so that V - E + F = 2
        
        self.num_nodes = 0
        self.num_edges = 0
        self.num_faces = 2
        
        # Create a dictionary whose keys are face labels, with sign to indicate
        # direction of rotation around face (CW = +, CCW = -); values are lists
        # of edges adjacent to face, and whose orientation is in given direction
        # around face
        
        self.face_dict = {}
        
        # Create a dictionary whose keys are face labels, with signs to indicate
        # rotation direction (as with self.face_dict); values are number of
        # edges in each list self.face_dict[face]
        
        self.edge_num_dict = {}
        
        # Create a dictionary whose keys are face labels *without* rotation
        # direction, and values are lists of edges in order going around the
        # face in CW order
        
        self.face_edge_dict = {}
        
        # Create a list of edge pair numbers used to reach this graph from the
        # primitive graph aa
        
        self.edge_pair_list = []
        
    def __repr__(self):
        
        symbol_list = []
        
        # See if there are any edges in the graph; if so, print them out as
        # ordered pairs.
        
        try:
            current = self.head.next
            
            while current.next:
                symbol_list.append(repr(current.first))
                second_symb = repr(current.last)
                
                current = current.next
                
            symbol_list.append(second_symb)
            
        except:
            return 'null'
            
        return ''.join(symbol_list)
    
    def initEdges(self):
        """
            This creates the graph with Gauss code aa in a standard form, and
            initializes the linked list by creating head, tail edges.
        """
        
        # Create symbols for node labels, and link together
        
        a_pos = Symbol(label = 1)
        a_neg = Symbol(label = -1)
        a_pos.linkSymbols(a_neg)
            
        # Standard notation for edges, and add head, tail to linked list
        
        edge_1 = Edge(first = a_pos, last = a_neg, face_left = 1, face_right = 2, \
                      edge_num = 0)
        edge_2 = Edge(first = a_neg, last = a_pos, face_left = 3, face_right = 1, \
                      prev = edge_1, next = None, edge_num = 1)
                      
        self.head = Edge(last = a_pos, prev = None, next = edge_1, edge_num = -1)
        
        edge_1.prev = self.head
        edge_1.next = edge_2
    
        # Initialize number of nodes, edges
        
        self.num_nodes = 1
        self.num_edges = 2
        self.num_faces = 3
        
        # Update face_dict and edge_num dictionaries
        
        for edge in [edge_1, edge_2]:
            self.face_dict[edge.face_left] = [edge]
            self.edge_num_dict[edge.face_left] = 1
            
            self.face_dict[-edge.face_right] = [edge]
            self.edge_num_dict[-edge.face_right] = 1
            
        # Update face_edge_dict
        
        self.face_edge_dict[1] = [edge_1, edge_2]
        self.face_edge_dict[2] = [edge_1]
        self.face_edge_dict[3] = [edge_2]
    
    def insertSelfLoop(self, edge_first = None, face_num = None):
        """
            This creates a self-loop inside the face face_num, by altering the
            edge edge_first, and adding two additional edges for the self-loop
            and the edge following it in the Eulerian circuit.
        """
            
        # Increment number of nodes and faces here, since node and face numbers
        # start from 1, not 0
        
        self.num_nodes += 1
        self.num_faces += 1
        
        # Get attributes of original edge, to use with other created edges
        
        last = edge_first.last
        next = edge_first.next
        
        face_left = edge_first.face_left
        face_right = edge_first.face_right
        
        index_left = self.face_edge_dict[face_left].index(edge_first)
        index_right = self.face_edge_dict[face_right].index(edge_first)
         
        # Create new final edge
            
        edge_last = Edge(last = last, face_left = face_left, face_right = face_right, \
                         next = next, edge_num = self.num_edges)
        
        if next:            # If an edge exists in linked list after edge_first
            next.prev = edge_last
        
        for face in [face_left, -face_right]:
            self.face_dict[face] = self.face_dict.get(face, []) + [edge_last]
            self.edge_num_dict[face] = self.edge_num_dict.get(face, 0) + 1
            
        # Create symbols for node labels, and link together
        
        a_pos = Symbol(label = self.num_nodes)
        a_neg = Symbol(label = -self.num_nodes)
        a_pos.linkSymbols(a_neg)
            
        # Place self-loop as appropriate, based on which side of original edge
        # the loop is placed
            
        if face_num == face_left:
            
            # Create self-loop
            
            edge_loop = Edge(first = a_neg, last = a_pos, prev = edge_first, \
                             next = edge_last, face_left = self.num_faces, \
                             face_right = face_num, edge_num = self.num_edges + 1)
            
            # Modify original edge by changing ending node; this node label is
            # opposite sign as starting node for edge_last
            
            edge_first.last = a_neg
            edge_last.first = a_pos
            
            # Update dictionaries
            
            for face in [self.num_faces, -face_left]:
                self.face_dict[face] = self.face_dict.get(face, []) + [edge_loop]
                self.edge_num_dict[face] = self.edge_num_dict.get(face, 0) + 1
                
            self.face_edge_dict[face_left] = self.face_edge_dict[face_left][:(index_left + 1)] + \
                [edge_loop, edge_last] + self.face_edge_dict[face_left][(index_left + 1):]
            self.face_edge_dict[face_right] = self.face_edge_dict[face_right][:index_right] + \
                [edge_last] + self.face_edge_dict[face_right][index_right:]
            
        elif face_num == -face_right:
            
            # Create self-loop
            
            edge_loop = Edge(first = a_pos, last = a_neg, prev = edge_first,
                             next = edge_last, face_left = -face_num, \
                             face_right = self.num_faces, edge_num = self.num_edges + 1)
            
            # Modify original edge by changing ending node; this node label is
            # opposite sign as starting node for edge_last
            
            edge_first.last = a_pos
            edge_last.first = a_neg
            
            # Update dictionaries
            
            for face in [face_right, -self.num_faces]:
                self.face_dict[face] = self.face_dict.get(face, []) + [edge_loop]
                self.edge_num_dict[face] = self.edge_num_dict.get(face, 0) + 1
                
            self.face_edge_dict[face_left] = self.face_edge_dict[face_left][:(index_left + 1)] + \
                [edge_last] + self.face_edge_dict[face_left][(index_left + 1):]
            self.face_edge_dict[face_right] = self.face_edge_dict[face_right][:index_right] + \
                [edge_last, edge_loop] + self.face_edge_dict[face_right][index_right:]
                
        # Add entry for self-loop in face_edge_dict
                
        self.face_edge_dict[self.num_faces] = [edge_loop]
            
        # Updating common to both scenarios: put self-loop between two edges,
        # and increase number of edges for edge_loop and edge_last
        
        edge_first.next = edge_loop
        edge_last.prev = edge_loop
        
        self.num_edges += 2
        
    def insertSymbol(self, edge_left = None, edge_right = None, face_num = None):
        
        # Increment number of nodes and faces here, since node and face numbers
        # start from 1, not 0
        
        self.num_nodes += 1
        self.num_faces += 1
        
        # Get symbols for each side of edges, and nodes on either side in linked
        # list
        
        last_left = edge_left.last
        next_left = edge_left.next
        
        first_right = edge_right.first
        prev_right = edge_right.prev
        
        # Get indices of edge_left, edge_right in face_edge_dict
        
        index_left = self.face_edge_dict[abs(face_num)].index(edge_left)
        index_right = self.face_edge_dict[abs(face_num)].index(edge_right)
        
        # Create new edge for "second half" of edge_left
        
        edge_begin = Edge(last = last_left, prev = edge_left, next = next_left, \
                          edge_num = self.num_edges)
            
        # Place edge_begin in linked list after edge_left, and before edge
        # originally after edge_left
            
        edge_left.next = edge_begin
        
        if next_left:   # If edge_left has an edge following it in linked list
            next_left.prev = edge_begin
        
        # Create new edge for "second half" of edge_right

        edge_end = Edge(first = first_right, prev = prev_right, next = edge_right, \
                        edge_num = self.num_edges + 1)
            
        # Place edge_end in linked list before edge_right, and after edge
        # originally before edge_right
        
        edge_right.prev = edge_end
        prev_right.next = edge_end
        
        # Create symbols for node labels, and link together
        
        a_pos = Symbol(label = self.num_nodes)
        a_neg = Symbol(label = -self.num_nodes)
        a_pos.linkSymbols(a_neg)
        
        # Whether the original edges have a CW or CCW rotation determines
        # which of the adjacent faces ends up on which side, as well as the
        # ordering of the signs for the node in the Gauss code. This info is
        # given by the sign of face_num.
        
        if face_num < 0:        # CW rotation (since common face is on right)
        
            # face_edge_dict only uses positive face numbers, and rotation of
            # edges already used, so switch sign
            
            face_num = -face_num
        
            # Put in appropriate node number for all edges; the last shall be
            # the first when begin and end are reversed, and we have left ->
            # end -> ... -> begin -> right
            
            edge_left.last = a_neg
            edge_end.last = a_neg
            
            edge_begin.first = a_pos
            edge_right.first = a_pos
        
            # Add in face information for begin, end edges
        
            edge_begin.face_left = edge_left.face_left
            edge_begin.face_right = self.num_faces
            
            edge_end.face_left = edge_right.face_left
            edge_end.face_right = self.num_faces
            
            self.face_dict[edge_left.face_left] += [edge_begin]
            self.face_dict[edge_right.face_left] += [edge_end]
            
            self.edge_num_dict[edge_left.face_left] += 1
            self.edge_num_dict[edge_right.face_left] += 1
            
            self.face_dict[-self.num_faces] = self.face_dict.get(-self.num_faces, []) + [edge_begin, edge_end]
            
            self.edge_num_dict[-self.num_faces] = self.edge_num_dict.get(-self.num_faces, 0) + 2
            
            # Update face_edge_dict
            
            cut_seq = [edge_begin, edge_end]
            
            if index_left < index_right:
                cut_seq += self.face_edge_dict[face_num][(index_right + 1):] + \
                    self.face_edge_dict[face_num][:index_left]
                    
                self.face_edge_dict[face_num] = self.face_edge_dict[face_num][index_left : (index_right + 1)]
            else:
                cut_seq += self.face_edge_dict[face_num][(index_right + 1) : index_left]
                
                self.face_edge_dict[face_num] = self.face_edge_dict[face_num][index_left:] + \
                    self.face_edge_dict[face_num][:(index_right + 1)]
                    
            # Put new edges in face_edge_dict before old edges; do this one at
            # a time, since the other faces may not be distinct
            
            index_begin = self.face_edge_dict[edge_left.face_left].index(edge_left)
            
            self.face_edge_dict[edge_left.face_left] = self.face_edge_dict[edge_left.face_left][:(index_begin + 1)] + \
                [edge_begin] + self.face_edge_dict[edge_left.face_left][(index_begin + 1):]
            
            index_end = self.face_edge_dict[edge_right.face_left].index(edge_right)
                
            self.face_edge_dict[edge_right.face_left] = self.face_edge_dict[edge_right.face_left][:index_end] + \
                [edge_end] + self.face_edge_dict[edge_right.face_left][index_end:]
            
        else:                   # CCW rotation (common face on left)
        
            # Put in appropriate node number for all edges; the last shall be
            # the first when begin and end are reversed, and we have left ->
            # end -> ... -> begin -> right
            
            edge_left.last = a_pos
            edge_end.last = a_pos
            
            edge_begin.first = a_neg
            edge_right.first = a_neg
        
            # Add in face information for begin, end edges
        
            edge_begin.face_left = self.num_faces
            edge_begin.face_right = edge_left.face_right
            
            edge_end.face_left = self.num_faces
            edge_end.face_right = edge_right.face_right
            
            self.face_dict[-edge_left.face_right] += [edge_begin]
            self.face_dict[-edge_right.face_right] += [edge_end]
            
            self.edge_num_dict[-edge_left.face_right] += 1
            self.edge_num_dict[-edge_right.face_right] += 1
            
            self.face_dict[self.num_faces] = self.face_dict.get(self.num_faces, []) + [edge_begin, edge_end]
            
            self.edge_num_dict[self.num_faces] = self.edge_num_dict.get(self.num_faces, 0) + 2
            
            # Update face_edge_dict
            
            cut_seq = [edge_end, edge_begin]
            
            if index_left < index_right:
                cut_seq += self.face_edge_dict[face_num][(index_left + 1) : index_right]
                
                self.face_edge_dict[face_num] = self.face_edge_dict[face_num][:(index_left + 1)] + \
                    self.face_edge_dict[face_num][index_right:]
            else:
                cut_seq += self.face_edge_dict[face_num][(index_left + 1):] + \
                    self.face_edge_dict[face_num][:index_right]
                    
                self.face_edge_dict[face_num] = self.face_edge_dict[face_num][index_right : (index_left + 1)]

            # Put new edges in face_edge_dict before old edges; do this one at
            # a time, since the other faces may not be distinct
            
            index_begin = self.face_edge_dict[edge_left.face_right].index(edge_left)
            
            self.face_edge_dict[edge_left.face_right] = self.face_edge_dict[edge_left.face_right][:index_begin] + \
                [edge_begin] + self.face_edge_dict[edge_left.face_right][index_begin:]
            
            index_end = self.face_edge_dict[edge_right.face_right].index(edge_right)
            
            self.face_edge_dict[edge_right.face_right] = self.face_edge_dict[edge_right.face_right][:(index_end + 1)] + \
                [edge_end] + self.face_edge_dict[edge_right.face_right][(index_end + 1):]

        # Create entry in face_edge_dict for new face, update adjacent faces for
        # affected edges in cut_seq
        
        self.face_edge_dict[self.num_faces] = cut_seq
        
        for edge in cut_seq:
            if edge.face_left == face_num:
                
                self.face_dict[edge.face_left].remove(edge)
                self.edge_num_dict[edge.face_left] -= 1
                
                edge.face_left = self.num_faces
                
                self.face_dict[self.num_faces] = self.face_dict.get(self.num_faces, []) + [edge]
                self.edge_num_dict[self.num_faces] = self.edge_num_dict.get(self.num_faces, 0) + 1
                
            elif edge.face_right == face_num:
                
                self.face_dict[-edge.face_right].remove(edge)
                self.edge_num_dict[-edge.face_right] -= 1
                
                edge.face_right = self.num_faces
                
                self.face_dict[-self.num_faces] = self.face_dict.get(-self.num_faces, []) + [edge]
                self.edge_num_dict[-self.num_faces] = self.edge_num_dict.get(-self.num_faces, 0) + 1
           
        # Reverse all edges between edge_left and edge_right
        
        self.reverse(left = edge_left, begin = edge_begin, end = edge_end, \
                      right = edge_right)
            
        # Clean up any entries in dictionaries with no elements; sort lists
        # with non-zero lengths. NOTE: A seperate list key_list is used here,
        # since entries may be deleted from edge_num_dict.keys
        
        key_list = [key for key in self.edge_num_dict.keys()]
        
        for key in key_list:
            if self.edge_num_dict[key] == 0:
                del self.edge_num_dict[key]
                del self.face_dict[key]
            else:
                self.face_dict[key].sort()
        
        # Increment number of edges
        
        self.num_edges += 2
        
        return [edge_left, edge_begin, edge_end, edge_right]
        
    def reverse(self, left = None, begin = None, end = None, right = None):
        """
            Reverse all symbols *between* left and right, starting at begin,
            and finishing at end. The symbols begin and end are kept separate
            in the reversal, since their node labels should already be correct
            from the generic splicing of the two original edges.
        """
        
        # Make sure that begin immediately follows left, and right immediately
        # follows end
        
        if left.next != begin or end.next != right:
            return False
        
        # See which edge comes first in linked list; if left comes before right,
        # flip order, so that reverse does not hit either end of linked list
        
        current = self.head.next
        while current != None:
            if current == left:
                break
            
            if current == right:
                left, begin, end, right = end, right, left, begin      
                break
            
            current = current.next
        
        # Shift ends of sequence in linked list
        
        left.next = end
        right.prev = begin
        
        # Reverse begin, but do not change sign for either symbol of the edge; these
        # should already be defined correctly from the step creating the new edges
        
        begin.prev = begin.next
        begin.next = right
        
        for face in [begin.face_left, -begin.face_right]:
            self.face_dict[face].remove(begin)
            self.edge_num_dict[face] -= 1
            
            self.face_dict[-face] = self.face_dict.get(-face, []) + [begin]
            self.edge_num_dict[-face] = self.edge_num_dict.get(-face, 0) + 1
        
        # Reverse signs of last node label
        
        begin.last.revSign()
        
        # Flip orientation
        
        begin.first, begin.last = begin.last, begin.first
        begin.face_left, begin.face_right = begin.face_right, begin.face_left
        
        current = begin.prev
        
        # Go through all edges between begin and end
        
        while current != end:
            
            prev_node = current.prev
           
            # Change placement of edge in face_dict
            
            for face in [current.face_left, -current.face_right]:
                self.face_dict[face].remove(current)
                self.edge_num_dict[face] -= 1
                
                self.face_dict[-face] = self.face_dict.get(-face, []) + [current]
                self.edge_num_dict[-face] = self.edge_num_dict.get(-face, 0) + 1
            
            # Reverse signs of last node label
            
            current.last.revSign()
            
            # Reverse order in linked list, and flip symbol order and face side
            
            current.prev = current.next
            current.next = prev_node
            
            current.first, current.last = current.last, current.first
            current.face_left, current.face_right = current.face_right, current.face_left
            
            # Move on to next edge
            
            current = current.prev
            
            # If the end of the linked list is reached, continue at beginning
            
            if current == None:
                current = self.head.next
        
        # Reverse end, but do not change sign for either symbol of the edge; these
        # should already be defined correctly from the step creating the new edges
        
        end.prev = left
        end.next = prev_node.prev
        
        for face in [end.face_left, -end.face_right]:
            self.face_dict[face].remove(end)
            self.edge_num_dict[face] -= 1
            
            self.face_dict[-face] = self.face_dict.get(-face, []) + [end]
            self.edge_num_dict[-face] = self.edge_num_dict.get(-face, 0) + 1
        
        end.first, end.last = end.last, end.first
        end.face_left, end.face_right = end.face_right, end.face_left
        
    def faceDict(self, face = None):
        if face:
            return repr(face) + ': ' + '[' + ', '.join([repr(edge) for edge in self.face_dict[face]]) + ']'

        return [repr(face) + ': ' + '[' + ', '.join([repr(edge) for edge in self.face_dict[face]]) + ']' \
            for face in sorted(self.face_dict.keys())]
                
    def faceEdgeDict(self, face = None):
        
        # If a specific face is asked for, return only the edges bordering that face
        
        if face:
            return repr(face) + ': ' + '[' + ', '.join([repr(edge) for edge in self.face_edge_dict[face]]) + ']'

        # For generic case, return all edges bordering all faces

        return [repr(face) + ': ' + '[' + ', '.join([repr(edge) for edge in self.face_edge_dict[face]]) + ']' \
            for face in sorted(self.face_edge_dict.keys())]
                
    def edgeNumDict(self, face = None):
        if face:
            return self.edge_num_dict.get(face, 0)

        return list(self.edge_num_dict.items())
        
    def edgeList(self):
        edge_list = []
        current = self.head.next
            
        while current != None:
            edge_list += [repr(current)]
            current = current.next
            
        return '[' + ', '.join(edge_list) + ']'
    
    def edgeEnds(self):
        current = self.head.next
        
        while current != None:
            yield repr(current.edge_num) + ': (' + repr(current.first) + ', ' + repr(current.last) + ')'
            current = current.next
        
    def numPairs(self, face = None):
        """
            Returns the total number of possible edge pairs for splicing. Edge
            pairs include (1) two distinct edges on the same face, and (2) two
            possible ways to splice an edge with itself, one for each face the
            edge is adjacent to. To avoid double-counting, these self-loops are
            counted once for each face. If the face is given, returns only the
            number of pairs for that face.
        """
        
        if face:
            edge_num = self.edge_num_dict.get(face, 0)
            return edge_num * (edge_num + 1) // 2
        
        return sum([iii * (iii + 1) // 2 for iii in self.edge_num_dict.values()])
        
    def getEdgePair(self, pair_num = None):
        """
            Return a pair of (not necessarily distinct) edges, and the face (with
            sign to indicate rotation around the face). If a number is given,
            return that number of pair of edges in list, if no number is specified,
            a random pair of edges is returned.
        """
        
        # Determine which face to get edges from
        
        tot_num_pairs = sum([iii * (iii + 1) // 2 for iii in self.edge_num_dict.values()])
        
        if type(pair_num) != int:
            pair_num = randrange(tot_num_pairs)
        elif pair_num >= tot_num_pairs:
            return False
        
        self.edge_pair_list += [pair_num]       # Add edge pair to list
        
        for face_num in sorted(self.edge_num_dict):
            num_edges = self.edge_num_dict[face_num]
            num_pairs = num_edges * (num_edges + 1) // 2
            pair_num -= num_pairs
            if pair_num < 0:
                break
            
        pair_num += num_pairs

        # Determine the edge(s) to choose
        
        iii = -1
        while pair_num >= 0:
            iii += 1
            jjj = pair_num
            
            pair_num -= (num_edges - iii)
            
        return [self.face_dict[face_num][iii], self.face_dict[face_num][iii + jjj], face_num]
    
    def allEdgePairs(self):
        """
            Creates an iterator which yields (edge_1, edge_2, face) for all
            possible edges on the same face that can be spliced together.
        """
        
        for face in self.edge_num_dict.keys():
            for pair in combinations_with_replacement(self.face_dict[face], 2):
                yield (pair[0], pair[1], face)
                
    def GaussToDT(self, f_list = False):
        """
            Converts the Gauss code into a DT sequence, as a list of node labels
            in the first entry of a list. If f_list is True, then the orientation
            information in the function f is returned as the second entry in the
            list. The resulting DT sequence is not guaranteed to be the
            lexicographically lowest form, so it is not necessarily in standard
            form.
        """
        
        DT_seq = [[] for iii in range(self.num_nodes)]
        f_list = [0 for iii in range(2 * self.num_nodes)]
        label = 0
        
        # Go through each edge, and convert Gauss code label to DT label
        
        current = self.head.next
        
        while current.next:
            second_symb = current.last
            
            if current.first.label > 0:
                DT_seq[current.first.label - 1] += [label]
                f_list[label] = +1
            else:
                DT_seq[-current.first.label - 1] += [label]
                f_list[label] = -1
                
            label += 1
            current = current.next
            
        if second_symb.label > 0:
            DT_seq[second_symb.label - 1] += [label]
            f_list[label] = +1
        else:
            DT_seq[-second_symb.label - 1] += [label]
            f_list[label] = -1
            
        # Sort each entry in DT_seq so that even number is first, then sort by
        # first label
        
        DT_seq = [[node[0], node[1]] if node[0] % 2 == 0 else [node[1], node[0]] for node in DT_seq]
        DT_seq.sort(key = lambda node : node[0])
        
        # Return results
            
        if f_list:
            return [DT_seq, f_list]
        else:
            return [DT_seq]
            
    def lowestCode(self):
        """
            For the graph, find the Gauss code with the lowest order. NOTE:
            This currently only works well for nine or less nodes, since in
            Python, '10...' < '2...'
        """
        
        symbol_list = []
        self_num_list = []
        
        try:
            current = self.head.next
            
            while current.next:
                
                self_num_list.append(current.first.label)
                second_num = current.last.label
                
                symbol_list.append(repr(current.first))
                second_symb = repr(current.last)
                
                current = current.next
                
            self_num_list.append(second_num)
            symbol_list.append(second_symb)
            
            # Create current string
            
            max_string = ''.join(symbol_list)
            
            for perm in permutations(range(1, self.num_nodes + 1)):
                
                # Try code in given order
                
                new_num_list = []
                
                for num in self_num_list:
                    if num > 0:
                        new_num_list += [str(perm[num - 1]) + '+']
                    else:
                        new_num_list += [str(perm[-num - 1]) + '-']
                        
                new_string = ''.join(new_num_list)
                
                if new_string < max_string:
                    max_string = new_string
                    
                # Try code in reversed order
                
                new_num_list = []
                
                for num in [self_num_list[iii] for iii in range(len(self_num_list) - 1, -1, -1)]:
                    if num > 0:
                        new_num_list += [str(perm[num - 1]) + '+']
                    else:
                        new_num_list += [str(perm[-num - 1]) + '-']
                        
                new_string = ''.join(new_num_list)
                
                if new_string < max_string:
                    max_string = new_string
            
        except:
            return ''
        
        return max_string
          
#-----------------------------------------------------------------------------#