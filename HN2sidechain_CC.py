#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 10:06:01 2024

@author: albertsmith
"""

import pyDR
import os
import numpy as np


md_dir='/Volumes/My Book/Y1/apo'

topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k}.xtc') for k in range(1,4)]

proj=pyDR.Project('Projects/HN2sidechain',create=True)

#%% iRED ana for chunks
from States import chunks

length=6600

for run,start,state,count in zip(*chunks(length).values()):
    sel=pyDR.MolSelect(topo,trajs[run],project=proj)
    sel.traj.t0=start
    sel.traj.tf=start+length
    sel.select_bond(Nuc='sidechain')
    i=[sel0.resnames[0]=='TRP' for sel0 in sel.sel1]
    sel1,sel2,rep=sel.sel1[i],sel.sel2[i],sel.repr_sel[i]
    sel.select_bond(Nuc='N')
    sel1=np.concatenate((sel1,sel.sel1))
    sel2=np.concatenate((sel2,sel.sel2))
    repr_sel=np.concatenate((rep,sel.repr_sel))
    sel.sel1,sel.sel2,sel.repr_sel=sel1,sel2,repr_sel
    
    pyDR.md2iRED(sel).iRED2data()
    proj[-1].source.additional_info=f'state{state}_{count}'
    
    
proj['iREDmode'].detect.r_no_opt(12)
proj['iREDmode'].fit().detect.r_auto(6)
proj['iREDmode']['no_opt'].fit().opt2dist(rhoz_cleanup=True)
proj['iREDmode']['opt_fit'].modes2bonds()

#%% Perform averaging
from pyDR.misc.Averaging import avgDataObjs
    
for state in range(1,4):
    avgDataObjs(proj[f'o6:IREDBOND:state{state}.'])
    
proj.save()
