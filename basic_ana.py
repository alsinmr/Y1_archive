#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 10:06:01 2024

@author: albertsmith
"""

import pyDR
import os


md_dir='/Volumes/My Book/Y1/apo'

topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k}.xtc') for k in range(1,4)]

proj=pyDR.Project('Projects/basic')


#%% Just HN for each run
for traj in trajs:
    sel=pyDR.MolSelect(topo,traj,project=proj)
    sel.traj.tf=15000
    sel.select_bond(Nuc='N')
    pyDR.md2data(sel)
    
proj.detect.r_no_opt(12)
proj.fit().detect.r_auto(6)
proj['no_opt'].fit().opt2dist(rhoz_cleanup=True)


#%% Chunks with just one state
length=6600  #Length of each ana
run=[1,1,2]  #Which sim?
start=[5900,12500,1170] #Starting frame
state=[2,2,2]