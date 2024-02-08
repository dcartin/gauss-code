# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 09:20:23 2019

@author: cartin
"""

def modOrbit(orig_perm, crossing = True, full_list = False):
    """
    Given a list of pairs of labels, where the first label is even and the
    second is odd, find the result of adding an integer b to all labels (mod
    the modulus) and adding it to the negative of the label. The procedure
    then returns all new non-isomorphic lists.
    
    The boolean variable crossing indicates whether the sequence has crossing
    information or not. If crossing is True, then it is assumed that the list
    is composed of three-item lists of the form (i, j, {-1, 0, 1}).
    
    When the odd and even labels are flipped, this also changes the signs of the
    crossings in the graph. (10 Nov 2019)
    
    Now only returns lowest lexicographic order graph, unless full_list is
    True. (05 Jul 2020)
    
    Compute modulus, rather than ask for it; I can't think of a reason not to
    simply use the obvious choice for all situations. Also made crossing = True
    the default, since that is usually what I use. (28 Sep 2023)
    """

    # # Assume that modulus is an even number

    # if modulus % 2 != 0:
    #     return False
    
    # Define modulus
    
    modulus = 2 * len(orig_perm)

    # We only need to look for all odd numbers, for example, less than the
    # value of the modulus, because then we can also add one to this for the
    # even numbers (not including the last one, which would be the value of the
    # modulus)

    bbb_list = range(1, modulus, 2)
    new_list = [[node for node in orig_perm]]
    
    # 27 Jul 2019: The algorithm below misses the transformation x -> -x. This
    # piece fixes that.

    if not crossing:
        temp_list = sorted([[(-pair[0]) % modulus, (-pair[1]) % modulus] for pair in orig_perm], key = lambda x : x[0])
    else:
        temp_list = sorted([[(-pair[0]) % modulus, (-pair[1]) % modulus, pair[2]] for pair in orig_perm], key = lambda x : x[0])

    if temp_list not in new_list:
        new_list += [temp_list]
            
    #-------------------------------------------------------------------------#
    
    if not crossing:
        temp_list = [sorted(sequence, key = lambda x : x[0]) for sequence in [[[(pair[1] + bbb) % modulus, (pair[0] + bbb) % modulus] for pair in orig_perm] for bbb in bbb_list]]
    else:
        temp_list = [sorted(sequence, key = lambda x : x[0]) for sequence in [[[(pair[1] + bbb) % modulus, (pair[0] + bbb) % modulus, -pair[2]] for pair in orig_perm] for bbb in bbb_list]]
    
    while len(temp_list) > 0:
        temp_seq = temp_list.pop()
        if temp_seq not in new_list:
            new_list += [temp_seq]
            
    #-------------------------------------------------------------------------#

    if not crossing:
        temp_list = [sorted(sequence, key = lambda x : x[0]) for sequence in [[[(bbb - pair[1]) % modulus, (bbb - pair[0]) % modulus] for pair in orig_perm] for bbb in bbb_list]]
    else:
        temp_list = [sorted(sequence, key = lambda x : x[0]) for sequence in [[[(bbb - pair[1]) % modulus, (bbb - pair[0]) % modulus, -pair[2]] for pair in orig_perm] for bbb in bbb_list]]

    while len(temp_list) > 0:
        temp_seq = temp_list.pop()
        if temp_seq not in new_list:
            new_list += [temp_seq]
            
    #-------------------------------------------------------------------------#

    # Now do the even numbers

    if not crossing:
        temp_list = [sorted(sequence, key = lambda x : x[0]) for sequence in [[[(pair[0] + bbb + 1) % modulus, (pair[1] + bbb + 1) % modulus] for pair in orig_perm] for bbb in bbb_list[:-1]]]
    else:
        temp_list = [sorted(sequence, key = lambda x : x[0]) for sequence in [[[(pair[0] + bbb + 1) % modulus, (pair[1] + bbb + 1) % modulus, pair[2]] for pair in orig_perm] for bbb in bbb_list[:-1]]]

    while len(temp_list) > 0:
        temp_seq = temp_list.pop()
        if temp_seq not in new_list:
            new_list += [temp_seq]
            
    #-------------------------------------------------------------------------#
    
    if not crossing:
        temp_list = [sorted(sequence, key = lambda x : x[0]) for sequence in [[[(bbb + 1 - pair[0]) % modulus, (bbb + 1 - pair[1]) % modulus] for pair in orig_perm] for bbb in bbb_list]]
    else:
        temp_list = [sorted(sequence, key = lambda x : x[0]) for sequence in [[[(bbb + 1 - pair[0]) % modulus, (bbb + 1 - pair[1]) % modulus, pair[2]] for pair in orig_perm] for bbb in bbb_list]]
        
    while len(temp_list) > 0:
        temp_seq = temp_list.pop()
        if temp_seq not in new_list:
            new_list += [temp_seq]
            
    #-------------------------------------------------------------------------#

    if full_list:
        return sorted(new_list)
    else:
        return sorted(new_list)[0]