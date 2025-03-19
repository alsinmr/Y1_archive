# -*- coding: utf-8 -*-

import pyDR
import os
from pyDR.PCA.PCAclean import PCA

import matplotlib.pyplot as plt
import numpy as np
from copy import copy
from pyDR.misc.tools import AA

from pyDR.Entropy import EntropyCC,CombineEntropy

TRP_res=[106,148,163,276,288]

step=1

proj=pyDR.Project('Projects/sidechain/')

md_dir='/Volumes/My Book/Y1/apo'

topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k+1}.xtc') for k in range(3)]

commands=['turn x -90','turn y 180','view','sel :276']

folder='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Report240603' #Location for  figures

select=pyDR.MolSelect(topo,trajs,project=proj)
select.select_bond('15N')
select.traj.step=step


ECC=EntropyCC(select)
ECC.load()


ECC.PCA=PCA(select)
ECC.PCA.select_atoms('name N C CA')
nc=6
ECC.PCA.Cluster.n_clusters=nc
ECC.PCA.Cluster.index=[0,1,2]
ECC.PCA.load()


# ax=ECC.plotCCpca()
# CCpca=ECC.CCpca
# a=np.argsort(CCpca)[-6:]
# for a0 in a:
#     ax.plot([a0,a0],CCpca[a0]+np.array([.01,.03]),color='red')
#     ax.text(a0,CCpca[a0]+.05,f'{AA(ECC.resi[a0].resname).symbol}{ECC.resids[a0]}',
#             rotation=90,horizontalalignment='center',verticalalignment='center')
# fig=ax.figure
# fig.set_size_inches([12.2,4.0])
# fig.tight_layout()


#%% Now loop over each separate trajectory
ECC=[ECC]
for traj in trajs:
    select=pyDR.MolSelect(topo,traj,project=proj)
    select.select_bond('15N')
    select.traj.step=step
    
    
    ECC.append(EntropyCC(select))
    ECC[-1].load()
    
    
    ECC[-1].PCA=PCA(select)
    ECC[-1].PCA.select_atoms('name N C CA')
    nc=6
    ECC[-1].PCA.Cluster.n_clusters=nc
    ECC[-1].PCA.Cluster.index=[0,1,2]
    ECC[-1].PCA.load()
    
# We'll just analyze the state of the first PCA, and replace PCA.Cluster.state with the results
l_apo=[] #Get the length of each run
for traj in trajs:
    select=pyDR.MolSelect(topo,traj)
    l_apo.append(len(select.traj))
    
for k,e in enumerate(ECC[1:]):
    e.PCA.Cluster.state=ECC[0].PCA.Cluster.state[sum(l_apo[:k]):sum(l_apo[:k+1])]
   
#  Plot all results
fig,ax=plt.subplots(2,1)
ECC[0].plotCCpca(ax=ax[0])

CCpca=ECC[0].CCpca
a=np.argsort(CCpca)[-6:]
for a0 in a:
    ax[0].plot([a0,a0],CCpca[a0]+np.array([.01,.03]),color='red')
    ax[0].text(a0,CCpca[a0]+.05,f'{AA(ECC[0].resi[a0].resname).symbol}{ECC[0].resids[a0]}',
            rotation=90,horizontalalignment='center',verticalalignment='center')
ax[0].set_ylim([0,.5])

for k,e in enumerate(ECC[1:]):
    e.plotCCpca(color=plt.get_cmap('tab10')(k),ax=ax[1])
ax[1].legend(('Run 1','Run 2','Run 3'))

fig.set_size_inches([7,7])
