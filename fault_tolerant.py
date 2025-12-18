#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 21:33:17 2023

@author: mark
"""

def gzerosgame(g,F=[],B=[]):
    """
    Return the derived set for a given graph g with set of banned edges B and a initial set of vertices. The derived set is given by doing generalized zero forcing process. That is, if y is the only white neighbor of x and xy is not banned, then x could force y into black.

    Input:
        g: a simple graph
        F: a list of vertices of g
        B: a list of tuples representing banned edges of g

    Output:
        A set of black vertices when zero forcing process stops.

    Examples:
        sage: gzerosgame(graphs.PathGraph(5),[0])
        set([0, 1, 2, 3, 4])
        sage: gzerosgame(graphs.PathGraph(5),[0],[(1,2)])
        set([0, 1])
    """
    S=set(F) # suspicuous vertices
    Black_vertices=set(F) # current black vertices
    again=1 # iterate again or not
    while again==1:
        again=0
        for x in S:
            N=set(g.neighbors(x))
            D=N.difference(Black_vertices) # set of white neighbors
            if len(D)==1:
                for v in D:
                    y=v # the only white neighbor
                if (((x,y) in B)==False) and (((y,x) in B)==False):
                    again=1
                    S.remove(x)
                    S.add(y)
                    Black_vertices.add(y)
                    break
    return(Black_vertices)

def ft_set(g,s, faults =1):
	# returns whether a set is a fault tolerant zero forcing set for a graph
    	# g a graph
	# s the set to test
	# faults the number of faults required
    Z=find_Z(g)
    if Z+faults>g.order():
        return False
    all_checks=True
    sub_s_count=0
    if len(zerosgame(g,s))==g.order():
            
        
        for sub_s in Subsets(s,len(s)-faults):
            
            #print("sub_s=",sub_s)
            
            if len(zerosgame(g,sub_s))==g.order():
                #print("it worked")
                sub_s_count+=1
                #print("subsets=",len(Subsets(s,len(s)-faults)))
                #print("count=",sub_s_count)
            if sub_s_count==len(Subsets(s,len(s)-faults)):
                return True    
            #else:
                #return False 
    return False

def ftZ(g,faults=1,robust=False,all_sets=False):
	#computes the fault tolerant zero forcing number
	#g a graph
	#faults the number of faults required
	#robust allows overfilling vertices, not implemented
	#all_sets=True returns all fault tolerant sets of minimum size
    ftz=-1
    Z=find_Z(g)
    ftz_sets=[]
    for i in range(faults,g.order()):
        
        S=[]
        for s in Subsets(g,Z+i):
            S.append(s)
            #print(S)
        for j in range(len(S)):
            if ft_set(g,S[j])==True:
                if ftz==-1 or len(S[j])==ftz:
                    ftz=len(S[j])
                    ftz_sets.append(S[j])
    if all_sets==True:
        return ftz_sets
    else:    
        return ftz	    

def psdgame(g,B=[],prop=False):
	#g a graph
	#blue a set of initially blue vertices
	again=1
	filled_vertices=set(B)
	vertices_to_add=set()
	num_rounds=0
	while again==1:
		again=0
		if len(filled_vertices)<g.order():
			num_rounds+=1
			#print(num_rounds)
			#print("filled_vertices =",filled_vertices )
		h=deepcopy(g)
		filled_vertices=filled_vertices.union(vertices_to_add)
		for b in filled_vertices:
			h.delete_vertex(b)
		hh=h.connected_components_subgraphs()
		for b in filled_vertices:
			for i in range(len(hh)):
				comp_vert=[]
				j=hh[i]
				comp_vert=[x for x in g.neighbors(b) if x in j.vertices() ]
				if len(comp_vert)==1:
					vertices_to_add.add(comp_vert[0])
					again=1
	if prop==True and len(filled_vertices)==g.order():
		return(num_rounds)
	else:
		return(filled_vertices)

def Zplus(g):
	#returns PSD forcing number of graph g
	size=1
	order=len(g.vertices())
	N=set(g.vertices())
	again=1
	while again==1:
		again=0
		for sub_s in Subsets(N,size):
			if len(psdgame(g,sub_s))==order:
				return size
		size+=1
		again=1    

def Z(g):
	size=1
	order=len(g.vertices())
	#print(order)
	N=set(g.vertices())
	again=1
	while again==1:
		again=0
		for sub_s in Subsets(N,size):
			#print(sub_s)
			#print(len(psdgame(g,sub_s)))
			if len(gzerosgame(g,sub_s))==order:
				#print("victory")
				return size
		#print("size=",size)
		size+=1
		again=1
def pt_plus(G,S):
    """
    pt_plus(G,S) returns the propagation time for
    positive semidefinite zero forcing on the graph G
    given a set S

    
    INPUT:
            graph G, vertex subset S
        
    OUTPUT: 
            Returns the propagation time based on the
            positive semidefinite color change rule.
    
    EXAMPLES:
        sage: pt_plus(graphs.CompleteGraph(6),[0,1,2,3])
        -1
        sage: pt_plus(graphs.PetersenGraph(),[0,2,8,7]) 
        2
    """
    prop_time = 0
    H = copy(G)
    H.delete_vertices(S)
    T = list(copy(S))
    while len(H.connected_components())==1:
        prop_time += 1
        new_T=copy(T)
        for v in T:
            N = set(G.neighbors(v))
            N.difference_update(set(T))
            if len(N) == 1:
                new_T.append(N.pop())
        if T == new_T:
            return -1
        T = new_T
        H = copy(G)
        H.delete_vertices(T)
    if len(H.connected_components())==0:
        return prop_time
    else:
        sub_times=[]
        for V in H.connected_components():
            B=set(G.vertices())
            B.difference_update(set(V))
            B=list(B)
            sub_times.append(pt_plus(G,B))
            if sub_times[-1]==-1:
                return -1
        return prop_time+max(sub_times)

def pt(G,S):
    """
    pt_plus(G,S) returns the propagation time for
    positive semidefinite zero forcing on the graph G
    given a set S

    
    INPUT:
            graph G, vertex subset S
        
    OUTPUT: 
            Returns the propagation time based on the
            positive semidefinite color change rule.
    
    EXAMPLES:
        sage: pt_plus(graphs.CompleteGraph(6),[0,1,2,3])
        -1
        sage: pt_plus(graphs.PetersenGraph(),[0,2,8,7]) 
        2
    """
    prop_time = 0
    H = copy(G)
    #H.delete_vertices(S)
    T = list(copy(S))
    while len(T)<G.order():
        prop_time += 1
        new_T=copy(T)
        #print("t=,",T)
        for v in T:
            N = set(G.neighbors(v))
            M=N.difference(set(T))
            #print("new_T=",new_T)
            if len(M) == 1:
                if list(M)[0] not in new_T:
                    new_T.append(M.pop())
                #else:
                    #print("error")
        #print(new_T)		
        if T == new_T:
            return -1
        T = new_T
    if len(T)==G.order():
        return prop_time

def twin(v,g):
    
    h=g.copy()
    nbrs=g.neighbors(v)
    w=g.order()
    h.add_vertex(w)
    for u in nbrs:
        h.add_edge(w,u)
    h.relabel(range(len(g.vertices()))) 
    return h

def is_minimal_zf_set(g,S):
    if len(zerosgame(g,S))<len(S):
        return False
    for s in Subsets(S,len(S)-1):
        if len(zerosgame(g,s))==g.order():
            return False
    return True        
        

def minimal_zf_sets(g):
    n=g.order()
    Z=find_Z(g)
    V=g.vertices()
    minzfsets=[]
    for i in range(Z,n):
        for subset in Subsets(V,i):
            if len(zerosgame(g,subset))==n and is_minimal_zf_set(g,subset):
                if s not in minzfsets:
                    minzfsets.append(subset)
    return minzfsets 

def minimum_zf_sets(g):
    n=g.order()
    Z=find_Z(g)
    V=g.vertices()
    minzfsets=[]
    
    for subset in Subsets(V,Z):
        if len(zerosgame(g,subset))==n and is_minimal_zf_set(g,subset):
            if subset not in minzfsets:
                minzfsets.append(subset)
    return minzfsets

def has_minimum_zf_set(g,S):
    newS=Set(S)
    zfsets=minimum_zf_sets(g)
    for s1 in zfsets:
        #print(type(s1))
        if Set(s1).issubset(newS):
            return True
    return False
# F-3 Function for propagation time of set S, pt(G,S)
# adapted by Leslie Hogben from Steve Butler's skew propagation time code (F-10)
# input: a graph G and a set S of vertices
# output: pt(G,S) = prop time of S in G (if -1 is returned then S is not a zero forcing set)
def ptz(G,S):
    V=set(G.vertices())
    count = 0
    done = False
    active = set(S)
    filled = set(S)
    for v in V:
        N=set(G.neighbors(v))
        if v in active and N.issubset(filled):
            active.remove(v)
        if (v in filled) and (v not in active) and (len(N.intersection(filled)) == G.degree(v)-1):
            active.add(v)   
    while not done:
        done = True
        new_active = copy(active)
        new_filled = copy(filled)
        for v in active:
            N=set(G.neighbors(v))
            if len(N.intersection(filled)) == G.degree(v)-1:
                if done:
                    done = False
                    count += 1
                N.symmetric_difference_update(N.intersection(filled))
                u=N.pop()
                new_active.remove(v)
                new_active.add(u)
                new_filled.add(u)
        active = copy(new_active)
        filled = copy(new_filled)
        # print filled
        for v in V:
            N=set(G.neighbors(v))
            if v in active and N.issubset(filled):
                active.remove(v)
            if (v in filled) and (v not in active) and (len(N.intersection(filled)) == G.degree(v)-1):
                active.add(v)
    if len(filled)==len(V):
        return count
    return -1

# F-4 Functions for standard propagation time 
# ptk(G,k) = min prop time over zero forcing sets of G of size k
# adapted by Leslie Hogben from Steve Butler's skew propagation time code (F-10)
# this is the function used to compute propation time of G and also throttling
# input: a graph G and a positive integer k >= Z(G)
# output: pt(G,k) = min prop time over zero forcing sets G of size k (if order+1 is returned then k<Z(G))
def ptk(G,k):
    V=set(G.vertices())
    sets = subsets(V,k)
    pt=len(V)+1
    for S in sets:
        count = 0
        done = False
        active = set(S)
        filled = set(S)
        for v in V:
            N=set(G.neighbors(v))
            if v in active and N.issubset(filled):
                active.remove(v)
            if (v in filled) and (v not in active) and (len(N.intersection(filled)) == G.degree(v)-1):
                active.add(v)   
        while not done:
            done = True
            new_active = copy(active)
            new_filled = copy(filled)
            for v in active:
                N=set(G.neighbors(v))
                if len(N.intersection(filled)) == G.degree(v)-1:
                    if done:
                        done = False
                        count += 1
                    N.symmetric_difference_update(N.intersection(filled))
                    u=N.pop()
                    new_active.remove(v)
                    new_active.add(u)
                    new_filled.add(u)
            active = copy(new_active)
            filled = copy(new_filled)
            # print filled
            for v in V:
                N=set(G.neighbors(v))
                if v in active and N.issubset(filled):
                    active.remove(v)
                if (v in filled) and (v not in active) and (len(N.intersection(filled)) == G.degree(v)-1):
                    active.add(v)
        if len(filled)==len(V):
            if count < pt:
                pt=count
    return pt
# propagation time pt(G)
# input: a graph G
# output: propagation time pt(G)=pt(G,Z(G))
def pt(g):
    ptg=ptk(g,Z(g))
    return ptg

# F-5 Functions for standard throttling number
#
# These functions compute th(G) by minimizing th(G,k)
#
# input: a graph G
# output: standard throttling number th(G) 
def th(g):
    z=Z(g)
    t=ptk(g,z)
    thz=z+t
    kmax=min(z+t-1,g.order())
    for k in [z+1..kmax]:
        ptkg=ptk(g,k)
        if k+ptkg<thz:
            thz=k+ptkg
    return thz
# input: a graph G
# output: the set [th(G),ko] such that th(G)= th(G,ko)
def thk(g):
    z=Z(g)
    ko=z
    t=ptk(g,z)
    thz=z+t
    kmax=min(z+t-1,g.order())
    for k in range(z,kmax+1):
        ptkg=ptk(g,k)
        if k+ptkg<thz:
            thz=k+ptkg
            ko=k
    return [thz,ko]

def ptz(G,S):
    V=set(G.vertices())
    count = 0
    done = False
    active = set(S)
    filled = set(S)
    for v in V:
        N=set(G.neighbors(v))
        if v in active and N.issubset(filled):
            active.remove(v)
        if (v in filled) and (v not in active) and (len(N.intersection(filled)) == G.degree(v)-1):
            active.add(v)   
    while not done:
        done = True
        new_active = copy(active)
        new_filled = copy(filled)
        for v in active:
            N=set(G.neighbors(v))
            if len(N.intersection(filled)) == G.degree(v)-1:
                if done:
                    done = False
                    count += 1
                N.symmetric_difference_update(N.intersection(filled))
                u=N.pop()
                new_active.remove(v)
                new_active.add(u)
                new_filled.add(u)
        active = copy(new_active)
        filled = copy(new_filled)
        # print filled
        for v in V:
            N=set(G.neighbors(v))
            if v in active and N.issubset(filled):
                active.remove(v)
            if (v in filled) and (v not in active) and (len(N.intersection(filled)) == G.degree(v)-1):
                active.add(v)
    if len(filled)==len(V):
        return count
    return -1

def ptk(G,k):
    V=set(G.vertices())
    sets = Subsets(V,k)
    pt=len(V)+1
    for S in sets:
        count = 0
        done = False
        active = set(S)
        filled = set(S)
        for v in V:
            N=set(G.neighbors(v))
            if v in active and N.issubset(filled):
                active.remove(v)
            if (v in filled) and (v not in active) and (len(N.intersection(filled)) == G.degree(v)-1):
                active.add(v)   
        while not done:
            done = True
            new_active = copy(active)
            new_filled = copy(filled)
            for v in active:
                N=set(G.neighbors(v))
                if len(N.intersection(filled)) == G.degree(v)-1:
                    if done:
                        done = False
                        count += 1
                    N.symmetric_difference_update(N.intersection(filled))
                    u=N.pop()
                    new_active.remove(v)
                    new_active.add(u)
                    new_filled.add(u)
            active = copy(new_active)
            filled = copy(new_filled)
            # print filled
            for v in V:
                N=set(G.neighbors(v))
                if v in active and N.issubset(filled):
                    active.remove(v)
                if (v in filled) and (v not in active) and (len(N.intersection(filled)) == G.degree(v)-1):
                    active.add(v)
        if len(filled)==len(V):
            if count < pt:
                pt=count
    return pt
# propagation time pt(G)
# input: a graph G
# output: propagation time pt(G)=pt(G,Z(G))
def pt(g):
    ptg=ptk(g,Z(g))
    return ptg

def ptpk(G,k):
    ord=G.order()
    V = G.vertices()
    S = Subsets(V,k)
    ptp = -1
    for s in S:
        ptps=pt_plus(G,s)
        if (ptp < 0):
            ptp=ptps
        if (ptps >= 0) and (ptps < ptp):
            ptp=ptps
    return ptp

def ptpk(G,S):
    ord=G.order()
    V = G.vertices()
    #S = Subsets(V,k)
    ptp = -1
    for s in S:
        ptps=pt_plus(G,s)
        if (ptp < 0):
            ptp=ptps
        if (ptps >= 0) and (ptps < ptp):
            ptp=ptps
    return ptp   

def ptpk(G,k):
    ord=G.order()
    V = G.vertices()
    S = Subsets(V,k)
    ptp = -1
    for s in S:
        ptps=pt_plus(G,s)
        if (ptp < 0):
            ptp=ptps
        if (ptps >= 0) and (ptps < ptp):
            ptp=ptps
    return ptp    
# function to compute pt_+(G) = pt_+(G,Z_+(G))
# input: a graph G  
# output: pt_+(G) 
def ptp(G):
    ptpg=ptpk(G,Zplus(G))
    return ptpg 

def is_minimal_pzf_set(g,S):
    if pt_plus(g,S)==-1:
        return False
    for s in Subsets(S,len(S)-1):
        if pt_plus(g,s)>=0:
            return False
    return True    

def minimal_pzf_sets(g):
    n=g.order()
    Z=Zplus(g)
    V=g.vertices()
    minpzfsets=[]

    for subset in Subsets(V,Z):
        if pt_plus(g,subset)>=0:
            for i in range(Z,n):
                for subset in Subsets(V,i):
                    if pt_plus(g,subset)>=0 and is_minimal_pzf_set(g,subset):
                        if subset not in minpzfsets:
                            minpzfsets.append(subset)
    return minpzfsets

def minimum_pzf_sets(g):
    n=g.order()
    Z=Zplus(g)
    V=g.vertices()
    minpzfsets=[]
    
    for subset in Subsets(V,Z):
        if pt_plus(g,subset)>=0 and is_minimal_pzf_set(g,subset):
            if subset not in minpzfsets:
                minpzfsets.append(subset)
    return minpzfsets    

def minimal_zf_size(g,interval=False):
    lengths=[]
    min_sets=minimal_zf_sets(g)
    for r in min_sets:
        if len(r) not in lengths:
            lengths.append(len(r))
    if interval==True:
        return lengths
    else:
        return max(lengths)

def has_fixed_prop_time(g):
    prop_times=[]
    zf_sets=minimal_zf_sets(g)
    for r in zf_sets:
        ptgs=ptz(g,r)
        if ptgs not in prop_times:
            prop_times.append(ptgs)
    if len(prop_times)==1 or 0:
        return True
    return False 

def minimal_pzf_size(g,interval=False):
    lengths=[]
    min_sets=minimal_pzf_sets(g)
    for r in min_sets:
        if len(r) not in lengths:
            lengths.append(len(r))
    if interval==True:
        return lengths
    else:
        return max(lengths)

def has_fixed_prop_time_plus(g):
    prop_times=[]
    zf_sets=minimal_pzf_sets(g)
    if len(zf_sets)==0:
        return True
    for r in zf_sets:
        ptgs=pt_plus(g,r)
        if ptgs not in prop_times:
            prop_times.append(ptgs)
    if len(prop_times)==1 or 0:
        return True
    return False    

def add_dom_vertex(g):
    h=g.copy()
    vert=h.order()
    for v in h.vertices():
        h.add_edge(vert,v)
    return h    

def threshold_graph(seq= [ 0]):
    n=len(seq)
    rng=range(n)
    g=graphs.EmptyGraph()
    for i in rng:
        verts=g.vertices()
        g.add_vertex(i)
        #g.show()
        if seq[i]!=0:
            for v in verts:
                g.add_edge(v,i)
    return g        

def LeakyForts(graph, leaks, PSD):
    """
    Returns a list of forts that can be used to determine leaky forcing number.
    
    A fort is appended to the list of forts if the total number of vertices outside the fort that are adjacent
    to exactly one vertex (or one vertex in each connected component, if PSD is True) 
    does not exceed the allowed number of leaks.
    
    Parameters:
      graph: A SageMath graph object.
      leaks: Number of leaks.
      PSD: if True, use PSD leaky forcing (process each connected component separately).
      
    Returns:
      A list of subsets (forts) that satisfy the leaky conditions, or a string message if too many leaks.
    """
    vertices = set(graph.vertices())
    if leaks >= len(vertices):
        return 'There are too many leaks'
    
    fort_list = []
    from sage.graphs.connectivity import connected_components_subgraphs
    
    # Iterate over all nonempty subsets of vertices as candidate forts.
    for fort in Subsets(list(vertices)):
        if not fort:
            continue
        non_fort = vertices - set(fort)
        leak_count = 0
        valid = True  # flag to exit early if leak_count exceeds leaks
        F = graph.subgraph(vertices=list(fort))
        
        if PSD:
            # Process each connected component of the fort
            for comp in connected_components_subgraphs(F):
                comp_vertices = set(comp.vertices())
                for u in non_fort:
                    count_deg = sum(1 for v in graph.neighbors(u) if v in comp_vertices)
                    if count_deg == 1:
                        leak_count += 1
                        if leak_count > leaks:
                            valid = False
                            break
                if not valid:
                    break
        else:
            # Standard leaky forcing: count neighbors in the entire fort
            fort_vertices = set(F.vertices())
            for u in non_fort:
                count_deg = sum(1 for v in graph.neighbors(u) if v in fort_vertices)
                if count_deg == 1:
                    leak_count += 1
                    if leak_count > leaks:
                        valid = False
                        break
        
        if valid and leak_count <= leaks:
            fort_list.append(fort)
    
    return fort_list

def LeakyNumber(graph, list_of_forts, leaks, Zflist, PSDforcing, k):
    """
    Computes the leaky forcing number for a graph given the a collection of forts.
    
    A set is a forcing set if it intersects every fort in list_of_forts.
    If Zflist is False, the function returns the size of a minimal forcing set (as an integer)
    as soon as one is found. If Zflist is True, it prints all forcing sets (though still returns
    the minimal forcing set size).
    
    Parameters:
      graph         : A SageMath graph.
      list_of_forts : A list of forts (each fort is an iterable of vertices).
      leaks         : Number of leaks.
      Zflist        : if True, prints the full list of forcing sets.
      PSDforcing    : if True, then it indicates whether PSD forcing is used.
      k             : If > 0, restricts search to candidate forcing sets of size exactly k.
                      If k is 0 or None, search proceeds by increasing subset size.
    
    Returns:
      The size (as an integer) of a minimal forcing set if found; otherwise, None.
    """
    # Check if the allowed leaks are too high compared to the maximum vertex degree.
    degrees = sorted(graph.degree())
    if leaks >= degrees[-1]:
        return(graph.order())  # All vertices would be in the leaky set.
    
    vertices = graph.vertices()
    
    # Pre-convert forts to sets for efficient membership testing.
    forts = [set(F) for F in list_of_forts]
    
    forcing_sets = []         # Will collect all forcing sets found.
    minimal_forcing_size = None
    found_min = False
    
    # Determine the sizes to iterate over:
    # If k > 0, consider only subsets of size k; otherwise, iterate from 1 to |V|.
    if k is not None and k > 0:
        subset_sizes = [k]
    else:
        subset_sizes = range(1, len(vertices) + 1)
    
    # Iterate over candidate forcing sets by increasing size.
    for r in subset_sizes:
        for candidate in Subsets(vertices, r):
            candidate_set = set(candidate)
            # A candidate is forcing if it intersects every fort.
            if all(candidate_set & F for F in forts):
                forcing_sets.append(candidate_set)
                if not found_min:
                    minimal_forcing_size = r
                    found_min = True
                    # If we don't need the full list, return immediately.
                    if not Zflist:
                        return int(minimal_forcing_size)
        # If we're not restricting by size and found a minimal set, no need to search larger subsets.
        if found_min and (k is None or k <= 0):
            break

    if not forcing_sets:
        return int(0) #return 0 for no forcing sets found

    if Zflist:
        print('The set of all ell-leaky sets is:', forcing_sets)
    
    return int(minimal_forcing_size)

