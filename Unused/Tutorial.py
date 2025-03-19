#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 09:48:52 2024

@author: albertsmith
"""

import pyDR
import os
import numpy as np
import matplotlib.pyplot as plt


#%% Part 0: Setup
pyDR.Defaults['zrange']=[-12,-1,200]

# File locations
md_dir='/Volumes/My Book/Y1/apo'
topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k+1}.xtc') for k in range(3)]

# Project (stores analyzed data)
proj=pyDR.Project()

# Selection: Specify bonds to examine
sel=pyDR.MolSelect(topo,trajs[0],project=proj)

sel.select_bond('N',fil)

# Set trajectory readout instructions (t0,tf, step)
sel.traj.step=10


#%% Part 1: Basic analysis
data=pyDR.md2data(sel,rank=1)

# Plot a correlation function
ax=plt.subplots()[1]
ax.plot(data.sens.info['t'],data.R[0,:])

data.detect.r_auto(6)
data.detect.plot_rhoz()

fit=data.fit()
opt=fit.opt2dist(rhoz_cleanup=True)

opt.plot(style='bar')
fit.plot()

#pyDR.chimeraX.chimeraX_funs.set_chimera_path('...') #Run this once

opt.chimera()

#%% Part 2: Cross-correlation analysis
ired=pyDR.md2iRED(sel,rank=1).iRED2data()

ired.detect.r_auto(6)
ifit=ired.fit()

iopt=ifit.opt2dist(rhoz_cleanup=True)


#%% Part 3: Frame analysis
helices=[np.arange(48,77),np.arange(84,113),np.arange(120,150),np.arange(165,185),
         np.arange(217,242),np.arange(261,292),np.arange(301,328),np.arange(332,343)]
for k in range(len(helices)):
    helices[k]-=6

sel.select_bond(Nuc='N',resids=np.concatenate(helices))
sel._mdmode=True
resids=sel.sel1.resids

for k,hlx in enumerate(helices):
    hlx0=[]
    for h in hlx:
        if h in resids:
            hlx0.append(h)
    helices[k]=np.array(hlx0)
    

frame_index=np.concatenate([np.ones(h.size)*k for k,h in enumerate(helices)])
sel0=[sel.sel1[np.isin(sel.sel1.resids,h)] for h in helices]

fo=pyDR.Frames.FrameObj(sel)
fo.rank=1

fo.tensor_frame(sel1=1,sel2=2)
fo.new_frame('peptide_plane',resids=sel.label)
fo.new_frame('superimpose',sel=sel0,frame_index=frame_index)
fo.new_frame('superimpose',sel=np.concatenate(sel0).sum(),frame_index=np.zeros(len(resids)))

fo.frames2data(mode='sym')

proj['Frames'].detect.r_auto(6)
proj['Frames'].fit().opt2dist(rhoz_cleanup=True)


#%% Part 4: PCA
sel.clear_sel()
pca=pyDR.PCA.PCA(sel)
pca.select_atoms('name CA N C')


