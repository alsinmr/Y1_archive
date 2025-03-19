#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 16:58:16 2024

@author: albertsmith
"""

import pyDR
import os
from pyDR.PCA.PCAclean import PCA

import matplotlib.pyplot as plt
import numpy as np
from copy import copy
from pyDR.misc.tools import AA

from pyDR.Entropy import EntropyCC,CombineEntropy

TRP_res=[106,148,163,276,288]

proj=pyDR.Project('Projects/sidechain/')

md_dir='/Volumes/My Book/Y1/apo'

topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k+1}.xtc') for k in range(3)]

commands=['turn x -90','turn y 180','view','sel :276']

folder='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/FigureParts' #Location for  figures

select=pyDR.MolSelect(topo,trajs,project=proj)
select.select_bond('15N')
select.traj.step=1


ECC=EntropyCC(select)
ECC.load()


ECC.PCA=PCA(select)
ECC.PCA.select_atoms('name N C CA')
nc=6
ECC.PCA.Cluster.n_clusters=nc
ECC.PCA.Cluster.index=[0,1,2]
fig=ECC.PCA.Cluster.plot(percent=False)[0].figure
fig.set_size_inches([8,3.6])
fig.tight_layout()

fig.savefig(os.path.join(folder,'PCA_cluster.png'))


#%% Chimera plots

ECC.CCchimera(indexCC='pca')

ECC.PCA.Hist.hist2struct()

#%% Plot the clusters
ECC.PCA.Cluster.plot()


def getW276state(i):
    from States import states
    states=np.concatenate([run for run in states.values()])
    return states[i]
    

W276state=getW276state(ECC.PCA.Cluster.cluster_index)

# Bar plot PCA state vs. Trp state
ax=plt.subplots()[1]
nt=ECC.PCA.Cluster.state.size
p=np.array([[(getW276state(ECC.PCA.Cluster.state==k)==q).sum()/nt for q in range(4)] for k in range(nc)])
hdl=[]
hatch=['','xxxx','////','oooo']
for k in range(4):
    for q in range(nc):
        color=[c for c in plt.get_cmap('tab10')(q)]
        # color[-1]=.25+k/4*.75
        h0=ax.bar(q,p[q,k:].sum()*100,color=color,hatch=hatch[k])
    hdl.append(h0)
ax.set_xlabel('PCA state')
ax.legend(hdl,('Unassigned','State 1','State 2','State 3'))
ax.set_title('W276 state vs. PCA state')
ax.set_ylabel('%')
fig=ax.figure
fig.set_size_inches([5.15,4.2])
fig.tight_layout()
fig.savefig(os.path.join(folder,'PCA_v_W276_percent.pdf'))

#%% Compare CC
ax=ECC.plotCCpca()
i=ECC.resids==276
ax.plot(ECC.CC[i].squeeze())
ax.set_xticks(np.arange(0,len(i),5))
