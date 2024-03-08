# -*- coding: utf-8 -*-


import numpy as np
import os
from MDAnalysis import Universe

states={f'run{k}':np.zeros(30000,dtype=int) for k in range(1,4)}


with open('/Volumes/My Book/Y1/apo/cheat_sheet_frame_chunks.txt') as f:
    for line in f:
        if np.any([f'run{k}' in line for k in range(1,4)]):
            run=line[:4]
            for k,frames in enumerate(line[4:].strip().split()):
                start,stop=[int(f) for f in frames.split('-')]
                states[run][start:stop+1]=k+1
            
md_dir='/Volumes/My Book/Y1/apo'
topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k}.xtc') for k in range(1,4)]

for k,traj in enumerate(trajs):
    uni=Universe(topo,traj)
    states[f'run{k+1}']=states[f'run{k+1}'][:len(uni.trajectory)]