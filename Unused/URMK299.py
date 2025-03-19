#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 09:52:35 2024

@author: albertsmith
"""

import os
from PCA import All as pca
import pyDR
from pyDR.Entropy import EntropyCC
import numpy as np
import matplotlib.pyplot as plt

fig_dir='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Report240828/figure_parts/anta'

proj=pyDR.Project()

#%% Ligand-Bound dynamics
from copy import copy

# Get clusters
pca2=pca.pcas[2]
pca2.Cluster.index=[2,3,4]
pca2.Cluster.n_clusters=5
pca2.Cluster.cluster_kwargs['random_state']=42
pca2.Cluster.plot()



# Time dependence of clusters
# Plot W276 state vs. time
cmap=plt.get_cmap('tab10')

fig,ax=plt.subplots(1,3)
l=pca2.traj.lengths
for k,a in enumerate(ax):
    state=pca2.Cluster.state[l[:k].sum():l[:k+1].sum()]
    
    a.plot(np.arange(state.size),state,color='black',linewidth=.5,linestyle=':')
    for m in range(pca2.Cluster.n_clusters):
        q=state==m
        a.scatter(np.arange(state.size)[q],np.ones(q.sum())*m,s=15,color=cmap(m))
        print(f'{q.sum()/q.__len__()*100:.1f}')
    a.set_xlabel('t / ns')
    a.set_yticklabels([])
    a.set_ylabel('state')
    a.set_xlim([-1300,17000])
    a.set_ylim([-.1,pca2.Cluster.n_clusters-0.9])
fig.set_size_inches([6.7*2,2.1*2])
fig.tight_layout()
fig.savefig(os.path.join(fig_dir,'anta_t_depend.png'))


# Cluster structures

fig=pca2.Cluster.plot(percent=False)[0].figure
ax=fig.axes
fig.set_size_inches([7.5,3.5])
fig.tight_layout()

proj.chimera.current=0
proj.chimera.close()
for k in range(pca2.Cluster.n_clusters):
    m=pca2.Cluster.cluster_index[k]
    pca2.traj[m]
    filename=os.path.join(fig_dir,f'anta_state{k+1}.pdb')
    pca2.uni.atoms.write(filename)
    clr=[int(c*100) for c in cmap(k)[:-1]]
    proj.chimera.command_line([f'open {filename}',f'color #{k+1} {clr[0]},{clr[1]},{clr[2]}',
                               'select ',f'contacts #{k+1}:276 restrict #{k+1}&~@CA,C,O,H*,N ignoreHiddenModels true color #000000'])
    if k:
        proj.chimera.command_line(f'align #{k+1}@CA toAtoms #1@CA')
        
    for q,a in enumerate(ax):
        a.scatter(pca2.PCamp[pca2.Cluster.index[q]][m],pca2.PCamp[pca2.Cluster.index[q+1]][m],color='black',marker='o')
        
fig.savefig(os.path.join(fig_dir,'anta_cluster.png'))


#%% Dynamics analysis
# Dynamic analysis of NPY bound

starts=[4000,10000,4000]
states=[1,2,5]
repeats=[0,0,0]
traj_index=[2,0,1]
length=6700

proj=pyDR.Project('Projects/URMK299_bound',create=True)

md_dir_as='/Volumes/My Book/Y1/Antagonist'

topo=os.path.join(md_dir_as,'Y1R_UR_MK299.pdb')
trajs=[os.path.join(md_dir_as,f'run{k+1}.xtc') for k in range(3)]

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
    
#%% 3D dynamics plots
hlx=[(36,68),(73,104),(109,144),(151,177),(204,246),(250,289),(295,323),(324,335)]

cmap=plt.get_cmap('gist_rainbow').resampled(10)
colors=[[c for c in cmap(k)] for k in range(9)]
nmdls=10
# 3D figures
proj.chimera.close()
for q,state in enumerate(['state1','state2','state5']):
    
    data=proj['iREDbond'][f'.+{state}'][0]
    
    frame=pca2.Cluster.cluster_index[int(state[-1])-1]
    if traj_index[q]==2:
        frame-=pca2.traj.lengths[:2].sum()
        
    elif traj_index[q]==1:
        frame-=pca2.traj.lengths[0]
    frame-=starts[q]
    frame=min(frame,len(data.select.traj)-1)
    data.select.molsys.make_pdb(ti=frame)
    
    data.select.molsys.chimera()
    
    index0=np.ones(data.label.shape,dtype=bool)
    i=np.argmax(data.label==276)
    for hlx0,color in zip(hlx,colors):
        resids=data.label
        index=np.logical_and(resids>=hlx0[0],resids<=hlx0[1])
        index0[index]=False
        data.select.chimera(index=index,x=data.CCnorm[-2][i][index],color=color)
    data.select.chimera(index=index0,x=data.CCnorm[-2][i][index0],color=colors[-1])
    
    if q:
        for m in range(nmdls):
            proj.chimera.command_line([f'align #{q*nmdls+m+1}@CA toAtoms #{m+1}@CA'])
        
proj.chimera.command_line(['~ribbon','~show','show :1-337','color #1 tan',
                           '~show ~@N,C,CA','turn x -90','view','turn y 135',
                           'color :276|:275@C,O,CA&~:276@C,O black','lighting full'])

for q,state in enumerate(['state1','state2','state5']):
    proj.chimera.command_line(['~show',f'show #{q*nmdls+1}-{(q+1)*nmdls}@N,C,CA',
                               f'show #{q*nmdls+1}-{(q+1)*nmdls}:276','~show @H*'])
    proj.chimera.savefig(os.path.join(fig_dir,f'anta276CC_{state}_rho4.png'),
                         options='transparentBackground True',overwrite=True)
    
    
hlx=[(36,68),(73,104),(109,144),(151,177),(204,246),(250,289),(295,323),(324,335)]

cmap=plt.get_cmap('gist_rainbow').resampled(10)
colors=[[c for c in cmap(k)] for k in range(10)]
# 3D figures
proj.chimera.close()
for q,state in enumerate(['state1','state2','state5']):
    
    data=proj['iREDbond'][f'.+{state}'][0]
    
    frame=pca2.Cluster.cluster_index[int(state[-1])-1]
    if traj_index[q]==2:
        frame-=pca2.traj.lengths[:2].sum()
        
    elif traj_index[q]==1:
        frame-=pca2.traj.lengths[0]
    frame-=starts[q]
    frame=min(frame,len(data.select.traj)-1)
    data.select.molsys.make_pdb(ti=frame)
    
    data.select.molsys.chimera()
    
    index0=np.ones(data.label.shape,dtype=bool)
    i=np.argmax(data.label==261)
    for hlx0,color in zip(hlx,colors):
        resids=data.label
        index=np.logical_and(resids>=hlx0[0],resids<=hlx0[1])
        index0[index]=False
        data.select.chimera(index=index,x=data.CCnorm[-2][i][index],color=color)
    data.select.chimera(index=index0,x=data.CCnorm[-2][i][index0],color=colors[-1])
    
    if q:
        for m in range(nmdls):
            proj.chimera.command_line([f'align #{q*nmdls+m+1}@CA toAtoms #{m+1}@CA'])
    
proj.chimera.command_line(['~ribbon','~show','show :1-337','color #1 tan',
                           '~show ~@N,C,CA','turn x -90','view','turn y 135',
                           'color :261|:260@C,O,CA&~:261@C,O black','lighting full'])
for q,state in enumerate(['state1','state2','state5']):
    proj.chimera.command_line(['~show',f'show #{q*nmdls+1}-{(q+1)*nmdls}@N,C,CA',
                               f'show #{q*nmdls+1}-{(q+1)*nmdls}:276','~show @H*'])    
    proj.chimera.savefig(os.path.join(fig_dir,f'anta261CC_{state}_rho4.png'),
                         options='transparentBackground True',overwrite=True)
    
    
    
hlx=[(36,68),(73,104),(109,144),(151,177),(208,241),(250,289),(295,323)]

# NPY amplitude (rho4)
cmap=plt.get_cmap('gist_earth').resampled(10)
colors=[[c for c in cmap(k)] for k in range(8)]
maxR=0
# 3D figures
proj.chimera.close()
nmdls=8
for q,state in enumerate(['state1','state2','state5']):
    
    data=proj['opt_fit']['MD'][f'.+{state}'][0]
    frame=pca2.Cluster.cluster_index[[1,0,3,2]][q]
    if traj_index[q]==2:
        frame-=pca2.traj.lengths[:2].sum()
        
    elif traj_index[q]==1:
        frame-=pca2.traj.lengths[0]
    frame-=starts[q]
    frame=min(frame,len(data.select.traj)-1)
    data.select.molsys.make_pdb(ti=frame)
    
    data.select.molsys.chimera()
    
    
    for hlx0,color in zip(hlx,colors):
        resids=data.label
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

for q,state in enumerate(['state1','state2','state5']):
    proj.chimera.command_line(['~show',f'show #{q*nmdls+1}-{(q+1)*nmdls}@N,C,CA'])
    proj.chimera.savefig(os.path.join(fig_dir,f'anta_rho4_{state}.png'),
                             options='transparentBackground True',overwrite=True)