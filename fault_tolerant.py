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
                print sub_s
                sub_s_count+=1
            else:
                all_checks=False
    return all_checks