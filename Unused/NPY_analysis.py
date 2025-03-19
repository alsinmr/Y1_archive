#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 14:04:40 2024

@author: albertsmith
"""

import os
from PCA import All as pca
import pyDR
from pyDR.Entropy import EntropyCC
import numpy as np
import matplotlib.pyplot as plt

fig_dir='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Report240828/figure_parts/npy'

proj=pyDR.Project()

#%% Ligand-Bound dynamics
from copy import copy

# Get clusters
pca1=pca.pcas[1]
pca1.Cluster.index=[0,1,2]
pca1.Cluster.n_clusters=4
pca1.Cluster.cluster_kwargs['random_state']=41
# Resort the states to be in order of activation
state=copy(pca1.Cluster.state)
state[pca1.Cluster.state==0]=2
state[pca1.Cluster.state==1]=0
state[pca1.Cluster.state==2]=3
state[pca1.Cluster.state==3]=1
pca1.Cluster._state=state


# Time dependence of clusters
# Plot W276 state vs. time
cmap=plt.get_cmap('tab10')

fig,ax=plt.subplots(1,3)
l=pca1.traj.lengths
for k,a in enumerate(ax):
    state=pca1.Cluster.state[l[:k].sum():l[:k+1].sum()]
    
    a.plot(np.arange(state.size),state,color='black',linewidth=.5,linestyle=':')
    for m in range(4):
        q=state==m
        a.scatter(np.arange(state.size)[q],np.ones(q.sum())*m,s=15,color=cmap(m))
        print(f'{q.sum()/q.__len__()*100:.1f}')
    a.set_xlabel('t / ns')
    a.set_yticklabels([])
    a.set_ylabel('state')
    a.set_xlim([-1300,17000])
    a.set_ylim([-.1,3.1])
fig.set_size_inches([6.7*2,2.1*2])
fig.tight_layout()
fig.savefig(os.path.join(fig_dir,'npy_state.png'))


md_dir_as='/Volumes/My Book/Y1/NPY'

topo=os.path.join(md_dir_as,'Y1R_NPY.pdb')
trajs=[os.path.join(md_dir_as,f'run{k+1}.proc.xtc') for k in range(3)]

commands=['turn x -90','turn y 180','view','sel :276']

select=pyDR.MolSelect(topo,trajs,project=proj)
select.select_bond('15N')
select.traj.step=1

ECCb=EntropyCC(select)

# Get structure for each cluster

fig=pca1.Cluster.plot(percent=False)[0].figure
ax=fig.axes
fig.set_size_inches([7.5,3.5])
fig.tight_layout()

proj.chimera.current=0
proj.chimera.close()
for k in range(4):
    m=pca1.Cluster.cluster_index[[1,3,0,2]][k]
    pca1.traj[m+100*(k==1)]
    filename=os.path.join(fig_dir,f'npy_state{k+1}.pdb')
    pca1.uni.atoms.write(filename)
    clr=[int(c*100) for c in cmap(k)[:-1]]
    proj.chimera.command_line([f'open {filename}',f'color #{k+1} {clr[0]},{clr[1]},{clr[2]}',
                               'select ',f'contacts #{k+1}:276 restrict #{k+1}&~@CA,C,O,H*,N ignoreHiddenModels true color #000000'])
    if k:
        proj.chimera.command_line(f'align #{k+1}@CA toAtoms #1@CA')
        
    for q,a in enumerate(ax):
        a.scatter(pca1.PCamp[q][m],pca1.PCamp[q+1][m],color='black',marker='o')
        
proj.chimera.command_line(['transparency 30 target rba','~show','show :276&~@H*'])

fig.savefig(os.path.join(fig_dir,'npy_cluster.png'))


# Dynamic analysis of NPY bound
starts=[6400,9000,0,10000]
states=[1,2,3,4]
repeats=[0,0,0,0]
traj_index=[1,0,2,2]
length=6700

proj=pyDR.Project('Projects/npy_bound',create=True)

md_dir_as='/Volumes/My Book/Y1/NPY'

topo=os.path.join(md_dir_as,'Y1R_NPY.pdb')
trajs=[os.path.join(md_dir_as,f'run{k+1}.proc.xtc') for k in range(3)]

if len(proj)<48:
    for start,traj,state,repeat in zip(starts,traj_index,states,repeats):
        select=pyDR.MolSelect(topo,trajs[traj],project=proj)
        select.traj.t0=start
        select.traj.tf=start+length
        select.select_bond('N')
        
        pyDR.md2data(select)
        pyDR.md2iRED(select,rank=1).iRED2data()
        
        title=f'BB_state{state}_rep{repeat}'
        proj[-2].source.additional_info=title
        proj[-1].source.additional_info=title
        
    proj.detect.r_no_opt(10)
    proj.fit().detect.r_auto(6)
    proj['no_opt'].fit().opt2dist(rhoz_cleanup=True).modes2bonds()
    
    proj.save()

# if len(proj)<52:
#     from pyDR.misc.Averaging import avgDataObjs
    
#     for sub in [proj['opt_fit']['MD'],proj['opt_fit']['iREDbond']]:
#         for state in ['state1','state2','state3','state4']:
#             avgDataObjs(sub[f'.+{state}'])
#             proj[-1].source.additional_info=proj[-1].source.additional_info.rsplit('_',maxsplit=1)[0]
        
#%% Plot dynamics

hlx=[(36,68),(73,104),(109,144),(151,177),(204,246),(250,289),(295,323),(324,335)]
hlx=[(36,68),(73,104),(109,144),(151,177),(208,241),(250,289),(295,323)]

# NPY amplitude (rho4)
cmap=plt.get_cmap('gist_earth').resampled(10)
colors=[[c for c in cmap(k)] for k in range(8)]
maxR=0
# 3D figures
proj.chimera.close()
nmdls=8
for q,state in enumerate(['state1','state2','state3','state4']):
    
    data=proj['opt_fit']['MD'][f'.+{state}'][0]
    data.select.molsys.chimera()
    frame=pca1.Cluster.cluster_index[[1,0,3,2]][q]
    if traj_index[q]==2:
        frame-=pca1.traj.lengths[:2].sum()
        
    elif traj_index[q]==1:
        frame-=pca1.traj.lengths[0]
    frame-=starts[q]
    data.select.traj[frame]
    
    for hlx0,color in zip(hlx,colors):
        resids=np.array([int(lbl.split('_')[1]) for lbl in data.label])
        index=np.logical_and(resids>=hlx0[0],resids<=hlx0[1])
        data.select.chimera(index=index,x=data.R[index][:,-2]*5,color=color)
        maxR=max(data.R[index][:,-2].max(),maxR)
        
    if q:
        for m in range(nmdls):
            proj.chimera.command_line([f'align #{q*nmdls+m+1}@CA toAtoms #{m+1}@CA'])
            
proj.chimera.command_line(['~ribbon','~show','show :1-337','color #1 tan',
                           '~show ~@N,C,CA,HN,O','turn x -90','view','turn y 150',
                           'show :276','color :276&~@N,C,CA,HN,O black',
                           'lighting full'])

for q,state in enumerate(['state1','state2','state3','state4']):
    proj.chimera.command_line(['~show',f'show #{q*nmdls+1}-{(q+1)*nmdls}@N,C,CA'])
    proj.chimera.savefig(os.path.join(fig_dir,f'npy_rho4_{state}.png'),
                             options='transparentBackground True',overwrite=True)
    
# NPY CC (rho4)
hlx=[(36,68),(73,104),(109,144),(151,177),(204,246),(250,289),(295,323),(324,335)]

cmap=plt.get_cmap('gist_rainbow').resampled(10)
colors=[[c for c in cmap(k)] for k in range(9)]
nmdls=11
# 3D figures
proj.chimera.close()
for q,state in enumerate(['state1','state2','state3','state4']):
    
    data=proj['iREDbond'][f'.+{state}'][0]
    data.select.molsys.chimera()
    frame=pca1.Cluster.cluster_index[[1,0,3,2]][q]
    if traj_index[q]==2:
        frame-=pca1.traj.lengths[:2].sum()
        
    elif traj_index[q]==1:
        frame-=pca1.traj.lengths[0]
    frame-=starts[q]
    data.select.traj[frame]
    
    index0=np.ones(data.label.shape,dtype=bool)
    i=np.argmax(data.label=='B_276')
    for hlx0,color in zip(hlx,colors):
        resids=np.array([int(lbl.split('_')[1]) for lbl in data.label])
        index=np.logical_and(resids>=hlx0[0],resids<=hlx0[1])
        index0[index]=False
        data.select.chimera(index=index,x=data.CCnorm[-2][i][index],color=color)
    index=np.array([lbl.split('_')[0]=='A' for lbl in data.label])
    index0[index]=False
    data.select.chimera(index=index,x=data.CCnorm[-2][i][index],color=color)
    data.select.chimera(index=index0,x=data.CCnorm[-2][i][index0],color=colors[-1])
    
    if q:
        for m in range(nmdls):
            proj.chimera.command_line([f'align #{q*nmdls+m+1}@CA toAtoms #{m+1}@CA'])
        
proj.chimera.command_line(['~ribbon','~show','show :1-337','color #1 tan',
                           '~show ~@N,C,CA','turn x -90','view','turn y 135',
                           'color :276|:275@C,O,CA&~:276@C,O black','lighting full'])

for q,state in enumerate(['state1','state2','state3','state4']):
    proj.chimera.command_line(['~show',f'show #{q*nmdls+1}-{(q+1)*nmdls}@N,C,CA',
                               f'show #{q*nmdls+1}-{(q+1)*nmdls}:276','~show @H*'])
    proj.chimera.savefig(os.path.join(fig_dir,f'npy276CC_{state}_rho4.png'),
                         options='transparentBackground True',overwrite=True)
    
    
hlx=[(36,68),(73,104),(109,144),(151,177),(204,246),(250,289),(295,323),(324,335)]

cmap=plt.get_cmap('gist_rainbow').resampled(10)
colors=[[c for c in cmap(k)] for k in range(10)]
# 3D figures
proj.chimera.close()
for q,state in enumerate(['state1','state2','state3','state4']):
    
    data=proj['iREDbond'][f'.+{state}'][0]
    data.select.molsys.chimera()
    
    frame=pca1.Cluster.cluster_index[[1,0,3,2]][q]
    if traj_index[q]==2:
        frame-=pca1.traj.lengths[:2].sum()
        
    elif traj_index[q]==1:
        frame-=pca1.traj.lengths[0]
    frame-=starts[q]
    data.select.traj[frame]
    
    index0=np.ones(data.label.shape,dtype=bool)
    i=np.argmax(data.label=='B_261')
    for hlx0,color in zip(hlx,colors):
        resids=np.array([int(lbl.split('_')[1]) for lbl in data.label])
        index=np.logical_and(resids>=hlx0[0],resids<=hlx0[1])
        index0[index]=False
        data.select.chimera(index=index,x=data.CCnorm[-2][i][index],color=color)
    index=np.array([lbl.split('_')[0]=='A' for lbl in data.label])
    index0[index]=False
    data.select.chimera(index=index,x=data.CCnorm[-2][i][index],color=color)
    data.select.chimera(index=index0,x=data.CCnorm[-2][i][index0],color=colors[-1])
    
    if q:
        for m in range(nmdls):
            proj.chimera.command_line([f'align #{q*nmdls+m+1}@CA toAtoms #{m+1}@CA'])
    
proj.chimera.command_line(['~ribbon','~show','show :1-337','color #1 tan',
                           '~show ~@N,C,CA','turn x -90','view','turn y 135',
                           'color :261|:260@C,O,CA&~:261@C,O black','lighting full'])
for q,state in enumerate(['state1','state2','state3','state4']):
    proj.chimera.command_line(['~show',f'show #{q*nmdls+1}-{(q+1)*nmdls}@N,C,CA',
                               f'show #{q*nmdls+1}-{(q+1)*nmdls}:276','~show @H*'])    
    proj.chimera.savefig(os.path.join(fig_dir,f'npy261CC_{state}_rho4.png'),
                         options='transparentBackground True',overwrite=True)
        