#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 21:33:17 2023

@author: mark
"""

def ft_set(g,s, faults =1):
    
    Z=find_gZ(g)
    if Z+faults>g.order():
        return -1
    all_checks=True
    if len(zerosgame(g,s))==g.order():
            
            
        for sub_s in subsets(s,len(s)-faults):
            sub_s_count=0
            if sub_s_count==len(sub_s):
                break
            if len(zerosgame(g,sub_s))==g.order():
                #print sub_s
                sub_s_count+=1
            else:
                all_checks=False
    return all_checks

def ftZ(g,faults=1,robust=False,all_sets=False):
    ftz=-1
    Z=find_gZ(g)
    ftz_sets=[]
    for i in range(faults,g.order()):
        S=[]
        for s in subsets(g,Z+i):
            S.append(s)
            
        for j in range(len(S)):
            if ft_set(g,S[j])==True:
                ftz=len(S[j])
                ftz_sets.append(S[j])
        if all_sets==True:
            return ftz_sets        
        else:    
            return ftz
