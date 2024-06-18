# -*- coding: utf-8 -*-

import pyDR
import os
from pyDR.PCA.PCAclean import PCA

import matplotlib.pyplot as plt
import numpy as np
from copy import copy
from pyDR.misc.tools import AA


step=1

TRP_res=[106,148,163,276,288]

proj=pyDR.Project('Projects/sidechain/')


# Load apo data
md_dir='/Volumes/My Book/Y1/apo'

topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k+1}.xtc') for k in range(3)]

l_apo=[]
for traj in trajs:
    select=pyDR.MolSelect(topo,traj)
    l_apo.append(len(select.traj))

commands=['turn x -90','turn y 180','view','sel :276']

folder='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Report240603' #Location for  figures

select=pyDR.MolSelect(topo,trajs,project=proj)
select.traj.step=step

apo=PCA(select)
apo.select_atoms('name N C CA and resid 28-339')
apo.align_ref='name CA and resid 28-339'
apo.load()



# Load bound data
#%% Overlap with activated state
md_dir_as='/Volumes/My Book/Y1/NPY'

topo=os.path.join(md_dir_as,'Y1R_NPY.pdb')
trajs=[os.path.join(md_dir_as,f'run{k+1}.proc.xtc') for k in range(3)]

l_bound=[]
for traj in trajs:
    select=pyDR.MolSelect(topo,traj)
    l_bound.append(len(select.traj))

select=pyDR.MolSelect(topo,trajs,project=proj)
select.traj.step=step

bound=PCA(select)
bound.select_atoms('name N C CA and segid B')

# Non-standard usage: we copy the reference positions from the apo trajectory
# and use them here so that we can align the two trajectories together.
bound._ref_pos = apo.ref_pos
bound.align_ref='name CA and segid B'
bound.load()

#%% Combine the two trajectories
both=copy(apo)
both._pos=np.concatenate((apo.pos,bound.pos),axis=0)

#Where are the two trajectories?
n0,n1=0,1
ax=both.Hist.plot(cmap='binary',n0=n0,n1=n1)

len_apo=apo.traj.__len__()

h=[]
for k in range(6):
    if k<3:
        i0=sum(l_apo[:k])
        i1=sum(l_apo[:k+1])
        color=plt.get_cmap('Blues')(.5+k*.25)
    else:
        i0=sum(l_apo)+sum(l_bound[:k-3])
        i1=sum(l_apo)+sum(l_bound[:k-2])
        color=plt.get_cmap('Reds')(.5+(k-3)*.25)
    h.append(ax.scatter(both.PCamp[n0][i0:i1],
                        both.PCamp[n1][i0:i1],s=1,color=color))

i0=27000
i1=33000
h.append(ax.scatter(both.PCamp[n0][i0:i1],
                    both.PCamp[n1][i0:i1],s=1,color='green'))
ax.legend(h,('apo1','apo2','apo3','bound1','bound2','bound3','closest approach'))
    

ax.legend(h,['apo','bound'])

# Where is the closest approach of the apo to the bound state?
boundPCmean=both.PCamp[:10][:,len_apo:].mean(1)
i=np.argmin(((both.PCamp[:10][:,:len_apo].T-boundPCmean)**2).sum(1))

i=np.argmax(both.PCamp[0][:len_apo])  #Yields 27778


def getW276state(i):
    from States import states
    states=np.concatenate([run for run in states.values()])
    return states[i]

getW276state(i)

# Trp state 2

nc=6
apo.Cluster.n_clusters=nc
apo.Cluster.index=[0,1,2]
apo.Cluster.state[i]