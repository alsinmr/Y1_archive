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
    select.traj.step=step
    l_apo.append(len(select.traj))

commands=['turn x -90','turn y 180','view','sel :276']

folder='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Report240603' #Location for  figures

select=pyDR.MolSelect(topo,trajs,project=proj)
select.traj.step=step

apo=PCA(select)
apo.select_atoms('name N C CA and resid 28-339')
apo.align_ref='name CA and resid 28-339'
apo.load()




#%% Load NPY-bound data
md_dir_as='/Volumes/My Book/Y1/NPY'

topo=os.path.join(md_dir_as,'Y1R_NPY.pdb')
trajs=[os.path.join(md_dir_as,f'run{k+1}.proc.xtc') for k in range(3)]

l_bound=[]
for traj in trajs:
    select=pyDR.MolSelect(topo,traj)
    select.traj.step=step
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

#%% Load antagonist-bound data
md_dir_as='/Volumes/My Book/Y1/Antagonist'

topo=os.path.join(md_dir_as,'Y1R_UR_MK299.pdb')
trajs=[os.path.join(md_dir_as,f'reduced_protein_1ns_center_helices_core_run{1}.xtc') for k in range(3)]

l_anta=[]
for traj in trajs:
    select=pyDR.MolSelect(topo,traj)
    select.traj.step=step
    l_anta.append(len(select.traj))

select=pyDR.MolSelect(topo,trajs,project=proj)
select.traj.step=step

anta=PCA(select)
anta.select_atoms('name N C CA and segid R and resid 28-339')

# Non-standard usage: we copy the reference positions from the apo trajectory
# and use them here so that we can align the two trajectories together.
anta._ref_pos = apo.ref_pos
anta.align_ref='name CA and segid R and resid 28-339'
anta.load()

#%% Load NPY/Gi-bound data
md_dir_as='/Volumes/My Book/Y1/NPY_Gi'

topo=os.path.join(md_dir_as,'Y1R_NPY_Gi.pdb')
trajs=[os.path.join(md_dir_as,f'run{1}.proc.xtc') for k in range(3)]

l_gi=[]
for traj in trajs:
    select=pyDR.MolSelect(topo,traj)
    select.traj.step=step
    l_gi.append(len(select.traj))

select=pyDR.MolSelect(topo,trajs,project=proj)
select.traj.step=step

gi=PCA(select)
gi.select_atoms('name N C CA and segid D and resid 28-339')

# Non-standard usage: we copy the reference positions from the apo trajectory
# and use them here so that we can align the two trajectories together.
gi._ref_pos = apo.ref_pos
gi.align_ref='name CA and segid D and resid 28-339'
gi.load()

#%% Combine the trajectories
    
from CombinedPCA import CombinePCA

All=CombinePCA(apo,bound,anta,gi)

three=CombinePCA(apo,bound,anta)

#Where are the two trajectories?
fig,ax=plt.subplots(2,2)
ax=ax.flatten()
for n0,n1,a in zip([0,1,2,3],[1,2,3,4],ax):
    # n0,n1=0,1
    All.Hist.plot(cmap='binary',n0=n0,n1=n1,ax=a)
    
    len_apo=apo.traj.__len__()
    
    h=[]
    for k in range(9):
        if k<3:
            i0=sum(l_apo[:k])
            i1=sum(l_apo[:k+1])
            color=plt.get_cmap('Blues')(.5+k*.25)
        elif k<6:
            i0=sum(l_apo)+sum(l_bound[:k-3])
            i1=sum(l_apo)+sum(l_bound[:k-2])
            color=plt.get_cmap('Reds')(.5+(k-3)*.25)
        elif k<9:
            i0=sum(l_apo)+sum(l_bound)+sum(l_anta[:k-6])
            i1=sum(l_apo)+sum(l_bound)+sum(l_anta[:k-5])
            color=plt.get_cmap('Greens')(.5+k*.25)
        elif k<12:
            i0=sum(l_apo)+sum(l_bound)+sum(l_anta)+sum(l_gi[:k-9])
            i1=sum(l_apo)+sum(l_bound)+sum(l_anta)+sum(l_gi[:k-8])
            color=plt.get_cmap('Purples')(.5+k*.25)
        h.append(a.scatter(All.PCamp[n0][i0:i1],
                            All.PCamp[n1][i0:i1],s=.1,color=color))

i0=27000//step
i1=33000//step
h.append(ax.scatter(All.PCamp[n0][i0:i1],
                    All.PCamp[n1][i0:i1],s=1,color='green'))
ax.legend(h,('apo1','apo2','apo3','bound1','bound2','bound3','anta1','anta2','anta3','gi1','gi2','gi3','closest approach'))
   

All.Cluster.index=[0,1,2,3,4,5]
All.Cluster.n_clusters=3
All.Cluster.plot()

# ax.legend(h,['apo','bound'])

# Where is the closest approach of the apo to the bound state?
boundPCmean=All.PCamp[:10][:,len_apo:].mean(1)
i=np.argmin(((All.PCamp[:10][:,:len_apo].T-boundPCmean)**2).sum(1))

i=np.argmax(All.PCamp[0][:len_apo])  #Yields 27778


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

#%% Some plotting
All.Cluster.index=[0,1,2,3,4]
All.Cluster.n_clusters=8
All.Cluster.plot()

if not(os.path.exists('pdbs')):os.mkdir('pdbs')

la=sum(l_apo)
lb=sum(l_bound)
for k,i in enumerate(All.Cluster.cluster_index):
    if i>la+lb:
        gi.traj[i-la-lb]
        gi.select.uni.atoms.select_atoms('segid D').write(f'pdbs/gi{i-la-lb}_cluster{k}.pdb')
    elif i>la:
        bound.traj[i-la]
        bound.select.uni.atoms.select_atoms('segid B').write(f'pdbs/bound{i-la}_cluster{k}.pdb')
    else:
        apo.traj[i]
        apo.select.uni.atoms.write(f'pdbs/apo{i}_cluster{k}.pdb')
        
        
proj.chimera.current=0
# proj.chimera.close()

k=0
for file in os.listdir('pdbs'):
    if '.pdb' in file:
        k+=1
        proj.chimera.command_line(f'open /Users/albertsmith/Documents/GitHub/Y1_archive/pdbs/{file}')
        color=[int(c*100) for c in plt.get_cmap('tab10')(int(file.split('cluster')[1][0]))[:3]]
        proj.chimera.command_line(f'color #{k} {color[0]},{color[1]},{color[2]}')
        if k>1:
            proj.chimera.command_line(f'align #{k}:28-339@CA toAtoms #1:28-339@CA')
        