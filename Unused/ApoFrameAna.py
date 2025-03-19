#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 14:16:32 2024

@author: albertsmith
"""


import pyDR
import os
import numpy as np


md_dir='/Volumes/My Book/Y1/apo'
topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k+1}.xtc') for k in range(3)]

proj=pyDR.Project('Projects/Frames',create=True)


select=pyDR.MolSelect(topo,trajs,project=proj)
select.traj.step=1


helices=[np.arange(48,77),np.arange(84,113),np.arange(120,150),np.arange(165,185),
         np.arange(217,242),np.arange(261,292),np.arange(301,328),np.arange(332,343)]
for k in range(len(helices)):
    helices[k]-=6

select.select_bond(Nuc='15N',resids=np.concatenate(helices))
select._mdmode=True
resids=select.sel1.resids

for k,hlx in enumerate(helices):
    hlx0=[]
    for h in hlx:
        if h in resids:
            hlx0.append(h)
    helices[k]=np.array(hlx0)
    

frame_index=np.concatenate([np.ones(h.size)*k for k,h in enumerate(helices)])
sel=[select.sel1[np.isin(select.sel1.resids,h)] for h in helices]

fo=pyDR.Frames.FrameObj(select)

fo.tensor_frame(sel1=1,sel2=2)
fo.new_frame('peptide_plane',resids=select.label)
fo.new_frame('superimpose',sel=sel,frame_index=frame_index)
fo.new_frame('superimpose',sel=np.concatenate(sel).sum(),frame_index=np.zeros(len(resids)))

fo.frames2data()

proj.detect.r_auto(6)
proj.fit().opt2dist(rhoz_cleanup=True)

