#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 21:33:17 2023

@author: mark
"""

def ft_set(g,s, faults =1):
    
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
                    print(S[j],"is fault tolerant")     
                    print(ftz_sets)                      
    if all_sets==True:
        return ftz_sets        
    else:    
        return ftz
