#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 16:02:37 2024

@author: albertsmith
"""

from pyDR.PCA.PCAclean import PCA
import pyDR
import os
import matplotlib.pyplot as plt
from copy import copy
import numpy as np


helices=[[37,68],[75,104],[110,144],[152,176],[206,244],[252,288],[296,323]]

md_dir='/Volumes/My Book/Y1/apo'

topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k}_whole.xtc') for k in range(1,4)]

step=1
sel=pyDR.MolSelect(topo,trajs,step=step)

#%% Backbone
BB=PCA(copy(sel))

BB.select_atoms('name N C CA')
BB.load()

#%% Helices only
helix=PCA(copy(sel))

helix.select_atoms('name N C CA and resid '+' '.join([f'{helix[0]}-{helix[1]}' for helix in helices]))
helix.load()

#%% Loops only
loop=PCA(copy(sel))

loop.select_atoms('name N C CA and not resid '+' '.join([f'{helix[0]}-{helix[1]}' for helix in helices]))
loop.load()

#%% Tryptophans
trp={}
for res in sel.uni.residues:
    if res.resname=='TRP':
        resid=str(res.resid)
        trp[resid]=PCA(copy(sel))
        trp[resid].select_atoms(f'resid {resid}')
        trp[resid].load()
    

#%% Plot comparison
fig,ax=plt.subplots(4,5)
for (key,value),ax0 in zip(trp.items(),ax.T):
    
    for pca,title,a in zip([BB,helix,loop,value],['Backbone','Helix only','Loops only',f'W{key}'],ax0):
        a.plot(pca.t,pca.PCamp[:2].T)
        a.set_title(title)
    ax0[-1].set_xlabel('t / ns')
fig.set_size_inches([15,8.8])
fig.tight_layout()


#%% Calculate CC for PC0
n0,n1=0,1
trp_CC={key:[] for key in trp}
for key,value in trp.items():
    for pca in [BB,helix,loop]:
        trp_CC[key].append((value.PCamp[n0]*pca.PCamp[n1]).mean()/np.sqrt(value.Lambda[n0]*pca.Lambda[n1]))
        
        
#%% Backbone with tryptophan (is the correlation?)
BBtrp=PCA(copy(sel))

BBtrp.select_atoms('name N C CA or resname TRP')
BBtrp.select_bond(Nuc='sidechain',filter_str='resname TRP')
BBtrp.load()
        