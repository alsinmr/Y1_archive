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

proj=pyDR.Project('Projects/basic',create=True)


#%% Just HN for each run
for traj in trajs:
    sel=pyDR.MolSelect(topo,traj,project=proj)
    sel.traj.tf=15000
    sel.select_bond(Nuc='N')
    pyDR.md2data(sel)
    



#%% Chunks with just one state
from States import chunks

length=6600
for run,start,state,count in zip(*chunks(length).values()):
    sel=pyDR.MolSelect(topo,trajs[run],project=proj)
    sel.traj.t0=start
    sel.traj.tf=start+length
    sel.select_bond(Nuc='N')
    pyDR.md2data(sel)
    proj[-1].source.additional_info=f'state{state}_{count}'
    
    
proj.detect.r_no_opt(12)
proj.fit().detect.r_auto(6)
proj['no_opt'].fit().opt2dist(rhoz_cleanup=True)

#%% iRED ana for chunks
length=6600
for run,start,state,count in zip(*chunks(length).values()):
    sel=pyDR.MolSelect(topo,trajs[run],project=proj)
    sel.traj.t0=start
    sel.traj.tf=start+length
    sel.select_bond(Nuc='N')
    pyDR.md2iRED(sel,rank=1).iRED2data()
    proj[-1].source.additional_info=f'state{state}_{count}'
    
    
proj['iREDmode'].detect.r_no_opt(12)
proj['iREDmode'].fit().detect.r_auto(6)
proj['iREDmode']['no_opt'].fit().opt2dist(rhoz_cleanup=True)
proj['iREDmode']['opt_fit'].modes2bonds()

#%% Perform averaging
from pyDR.misc.Averaging import avgDataObjs

for state in range(1,4):
    avgDataObjs(proj[f'o6:MD:state{state}.'])
    
for state in range(1,4):
    avgDataObjs(proj[f'o6:IREDBOND:state{state}.'])
    
proj.save()