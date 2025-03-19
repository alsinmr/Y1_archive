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
    

class StateChunk():
    def __init__(self,states):
        self.states=states
        
    def __call__(self,length=6000):
        states=list()
        runs=list()
        starts=list()
        counts=list()
        
        for k,x in enumerate(self.states.values()):
            pos=0
            while pos<len(x):
                if pos+length<=len(x) and np.all(x[pos]==x[pos:pos+length]):
                    runs.append(k)
                    starts.append(pos)
                    states.append(x[pos])
                    counts.append(np.sum([s==states[-1] for s in states[:-1]]).astype(int))
                    pos+=length
                else:
                    q=np.diff(x[pos:])!=0
                    if np.any(q):
                        pos+=np.argmax(q)+1
                    else:
                        pos=len(x)
        
        out={'run':runs,'start':starts,'state':states,'count':counts}
        
        i=np.argsort(states)
        out={key:np.array(value)[i] for key,value in out.items()}
        i=out['state']!=0
        out={key:value[i] for key,value in out.items()}
        
        return out
    
chunks=StateChunk(states)
                    
        