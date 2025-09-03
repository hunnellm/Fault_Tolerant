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
