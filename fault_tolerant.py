#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 21:33:17 2023

@author: mark
"""

def ft_set(g,s, faults =1):
	# returns whether a set is a fault tolerant zero forcing set for a graph
    	# g a graph
	# s the set to test
	# faults the number of faults required
    Z=find_Z(g)
    if Z+faults>g.order():
        return -1
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

def psdgame(g,B=[]):
	#returns the set of vertices colored by the psd color change tule
	#g a graph
	#B the initially blue vertices
	again=1
	filled_vertices=set(B)
	vertices_to_add=set()
	while again==1:
		again=0
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
    else:    
        return ftz
