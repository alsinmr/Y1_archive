#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 09:03:50 2024

@author: albertsmith
"""

import os
from ClusterAna import ECC,pca
import pyDR
from pyDR.Entropy import EntropyCC
import numpy as np
import matplotlib.pyplot as plt


# Figure directory
fig_dir='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Report240828/figure_parts/apo'



#%% pca clustering is performed in ClusterAna
proj=pca.project
    

# Backbone clustering    
figBB=pca.pcas[0].Cluster.plot(percent=False)[0].figure
figBB.set_size_inches([7.5,7])
figBB.tight_layout()
figBB.savefig(os.path.join(fig_dir,'apo_cluster.png'))

# Correlation to sidechains

index=np.argmax(ECC.resids==276)

N=ECC.index[3:,index].sum()

chi=[ECC.chi[k][ECC.index[k+3,:index].sum()] for k in range(N-1,-1,-1)]
from sklearn import cluster
from sklearn.cluster import MiniBatchKMeans
cluster=MiniBatchKMeans(n_clusters=4,random_state=3)
state=cluster.fit(np.array(chi).T).labels_
ECC.state[index]=state
ECC.plotChi(index)

ax=ECC.plotCCpca()
fig=ax.figure
ax.set_xticks(np.arange(0,len(ECC.resi),10))
fig.set_size_inches([12.5,4])
fig.tight_layout()

fig.savefig(os.path.join(fig_dir,'sidechain_pca_CC.pdf'))

# Plot sidechain entropy
fig,ax=plt.subplots()
ax.plot(np.arange(len(ECC.resids)),ECC.Sres)
ax.xaxis.set_major_formatter(ECC._axis_formatter(index=np.ones(ECC.resids.shape,dtype=bool)))

ax.set_xticks(np.arange(0,len(ECC.resids),10))
ax.tick_params(axis='x', labelrotation=90)
ax.set_ylabel(r'S$_{res}$ / J*mol$^{-1}$*K$^{-1}$')
fig.set_size_inches([12.5,4])
fig.tight_layout()
fig.savefig(os.path.join(fig_dir,'sidechain_entropy.pdf'))

# PCA vs. 276W state
pc_state=[(86,-86,0),(66,-48,1),(30,-23,2),(55,2,3)]

fig,ax=plt.subplots(1,2)
ax=ax.flatten()
step=1
cmap=plt.get_cmap('turbo').resampled(5)
i=np.argmax(ECC.resids==276)

# proj.chimera.current=0
# proj.chimera.close()

pc_state_ind=[]


for n0,n1,a in zip(np.arange(0,len(ax)),np.arange(1,len(ax)+1),ax):
    pca.pcas[0].Hist.plot(n0,n1,cmap='Reds',ax=a)
    for k in range(4):
        q=ECC.state[i]==k
        print(f'State {k}: {q.sum()/len(q)*100:.1f} %')
        a.scatter(pca.pcas[0].PCamp[n0][q][::step],pca.pcas[0].PCamp[n1][q][::step],color=cmap(k+1),s=.1)
    for n,(pc0,pc1,state) in enumerate(pc_state):
        q=ECC.state[i]==state
        m=np.argmin((pca.pcas[0].PCamp[0][q]-pc0)**2+(pca.pcas[0].PCamp[1][q]-pc1)**2)
        m=np.argwhere(q)[:,0][m]
        if n0==0:pc_state_ind.append(m)
        x,y=pca.pcas[0].PCamp[[n0,n1]][:,m]
        a.scatter(x,y,color=cmap(n+1),marker=['o','^','*','s'][n],s=40,edgecolors='black')
    a.set_xlim([-130,130])
    a.set_ylim([-130,130])
        
fig.set_size_inches([6.8,3.2])
fig.tight_layout()
fig.savefig(os.path.join(fig_dir,'apoPCA_vs_W276state.png'))



fig,ax=plt.subplots(1,2)
step=1
i=np.argmax(ECC.resids==276)
cmaps=['Blues','Greens','Oranges','Greys']

from matplotlib.colors import ListedColormap

for n0,n1,a in zip([0,1],[1,2],ax):
    # pca.pcas[0].Hist.plot(n0,n1,cmap='Reds',ax=a)
    for n,((*_,k),cmap) in enumerate(zip(pc_state,cmaps)):
        colors=[[*[c for c in plt.get_cmap(cmap)(k)[:3]],(k-100)/156] for k in range(100,256)]
        colors[0]=(0,0,0,0)
        cmap=ListedColormap(colors)
        q=ECC.state[i]==k
        print(f'State {k}: {q.sum()/len(q)*100:.1f} %')
        # a.scatter(pca.pcas[0].PCamp[n0][q][::step],pca.pcas[0].PCamp[n1][q][::step],color=cmap(n+1),s=.1)
        a.hist2d(pca.pcas[0].PCamp[n0][q][::step],pca.pcas[0].PCamp[n1][q][::step],cmap=cmap,bins=50)
    for n,(pc0,pc1,state) in enumerate(pc_state):
        q=ECC.state[i]==state
        m=np.argmin((pca.pcas[0].PCamp[0][q]-pc0)**2+(pca.pcas[0].PCamp[1][q]-pc1)**2)
        m=np.argwhere(q)[:,0][m]
        if n0==0:pc_state_ind.append(m)
        x,y=pca.pcas[0].PCamp[[n0,n1]][:,m]
        a.scatter(x,y,color=cmap(n+1),marker=['o','^','*','s'][n],s=40,edgecolors='black')
    
    a.set_xlim([-130,130])
    a.set_ylim([-130,130])
    a.set_xlabel(f'PC {n0}')
    a.set_ylabel(f'PC {n1}')

fig.set_size_inches([6.8,3.2])
fig.tight_layout()

# Add scatter points to backbone clustering plot
markers=['*','^','p','P']
for k,a in enumerate(figBB.axes):
    for q,i in enumerate(pc_state_ind[:-1]):
        a.scatter(pca.pcas[0].PCamp[k][i],pca.pcas[0].PCamp[k+1][i],s=50,
                  color=cmap(q+1),marker=markers[q],edgecolor='black')
figBB.savefig(os.path.join(fig_dir,'apo_cluster.png'))

# Plot states in chimeraX

# PCA vs. W276 bar plot
ax=plt.subplots()[1]
hatch=['xxx','///','ooo',''][::-1]
for k in range(6):
    count0=pca.pcas[0].Cluster.state==k
    count1=np.zeros(count0.shape,dtype=bool)
    y=[]
    for m,(*_,state) in enumerate(pc_state):
        count1=np.logical_or(count1,ECC.state[i]==state)
        y.append(np.logical_and(count0,count1).sum()/len(count1)*100)
    for m,y0 in enumerate(y[::-1]):
        ax.bar(k,y0,color=plt.get_cmap('tab10')(k),hatch=hatch[m])
ax.legend(('state 4','state 3','state 2','state 1'))
fig=ax.figure
fig.set_size_inches([7.5,7])
fig.savefig(os.path.join(fig_dir,'PCAcluster_v_W276state.png'))

# Get structures
# We want to capture both certain areas in the PCA, but also certain tryptophan
# states

proj.chimera.current=0
proj.chimera.close()
for k,(pc0,pc1,state) in enumerate(pc_state):
    q=ECC.state[i]==state
    m=np.argmin((pca.pcas[0].PCamp[0][q]-pc0)**2+(pca.pcas[0].PCamp[1][q]-pc1)**2)
    m=np.argwhere(q)[:,0][m]
    pca.pcas[0].traj[m]
    filename=os.path.join(fig_dir,f'Trp_state{state}.pdb')
    pca.pcas[0].uni.atoms.write(filename)
    clr=[int(c*100) for c in cmap(k+1)[:-1]]
    proj.chimera.command_line([f'open {filename}',f'color #{k+1} {clr[0]},{clr[1]},{clr[2]}',
                               'select ',f'contacts #{k+1}:276 restrict #{k+1}&~@CA,C,O,H*,N ignoreHiddenModels true color #000000'])
    if k:
        proj.chimera.command_line(f'align #{k+1}@CA toAtoms #1@CA')

proj.chimera.command_line(['transparency 30 target rba','~show','show :276'])

# Plot CC onto molecule
ECC.CCchimera(indexCC='pca')
    
    
# Plot W276 state vs. time
fig,ax=plt.subplots(1,3)
l=ECC.traj.lengths
for k,a in enumerate(ax):
    state0=ECC.state[i,l[:k].sum():l[:k+1].sum()]
    state=np.zeros(state0.shape)
    state[state0==4]=1
    state[state0==0]=2
    state[state0==1]=3
    state[state0==5]=4
    
    a.plot(np.arange(state.size),state,color='black',linewidth=.5,linestyle=':')
    for m in range(1,5):
        q=state==m
        a.scatter(np.arange(state.size)[q],np.ones(q.sum())*m,s=15,color=cmap(m))
        print(f'{q.sum()/q.__len__()*100:.1f}')
    a.set_xlabel('t / ns')
    a.set_yticklabels([])
    a.set_ylabel('state')
    a.set_xlim([-1300,27600])
fig.set_size_inches([8.7,2.1])
fig.tight_layout()
fig.savefig(os.path.join(fig_dir,'Trp276_state.png'))


#%% Apo dynamics
md_dir='/Volumes/My Book/Y1/apo'
topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k+1}.xtc') for k in range(3)]

proj=pyDR.Project('Projects/apo',create=True)
states=[1,1,1,2,2,2,3]
repeats=[0,1,2,0,1,2,0]
traj_index=[1,1,2,0,0,1,0]
starts=[14000,9567,0,6000,12700,1250,1180,19160]
length=6700

if len(proj)<49:
    for start,traj,state,repeat in zip(starts,traj_index,states,repeats):
        select=pyDR.MolSelect(topo,trajs[traj],project=proj)
        select.traj.t0=start
        select.traj.tf=start+length
        select.select_bond('N')
        
        pyDR.md2data(select)
        pyDR.md2iRED(select,rank=1).iRED2data()
        
        title=f'W276_state{state}_rep{repeat}'
        proj[-2].source.additional_info=title
        proj[-1].source.additional_info=title
        
    proj.detect.r_no_opt(10)
    proj.fit().detect.r_auto(6)
    proj['no_opt'].fit().opt2dist(rhoz_cleanup=True).modes2bonds()
    
    proj.save()

if len(proj)<55:
    from pyDR.misc.Averaging import avgDataObjs
    
    for sub in [proj['opt_fit']['MD'],proj['opt_fit']['iREDbond']]:
        for state in ['state1','state2','state3']:
            avgDataObjs(sub[f'.+{state}'])
            proj[-1].source.additional_info=proj[-1].source.additional_info.rsplit('_',maxsplit=1)[0]
        

hlx=[(36,68),(73,104),(109,144),(151,177),(204,246),(250,289),(295,323),(324,335)]

hlx=[(36,68),(74,102),(110,144),(154,175),(209,241),(257,289),(295,323)]
index=np.zeros(proj[-1].label.shape,dtype=bool)
for hlx0 in hlx:
    i=np.logical_and(proj[-1].label>=hlx0[0],proj[-1].label<=hlx0[1])
    index[i]=True
    
# proj['MD']['.+AvOb'].chimera(scaling=5,index=index)


# Apo amplitude (rho4)
cmap=plt.get_cmap('gist_earth').resampled(10)
colors=[[c for c in cmap(k)] for k in range(8)]
maxR=0
# 3D figures
for state in ['state1','state2','state3']:
    proj.chimera.close()
    data=proj['MD'][f'.+AvOb.+{state}'][0]
    data.select.molsys.chimera()
    for hlx0,color in zip(hlx,colors):
        index=np.logical_and(data.label>=hlx0[0],data.label<=hlx0[1])
        data.select.chimera(index=index,x=data.R[index][:,-2]*10,color=color)
        maxR=max(data.R[index][:,-2].max(),maxR)
        
    proj.chimera.command_line(['~ribbon','~show','show :1-337','color #1 tan',
                               '~show ~@N,C,CA,HN,O','turn x -90','view','turn y 135',
                               'show :276','color :276&~@N,C,CA,HN,O black',
                               'lighting full'])
    
    proj.chimera.savefig(os.path.join(fig_dir,f'apo_{state}_rho4_new.png'),
                         options='transparentBackground True',overwrite=True)
    
# Apo CC (rho4)
hlx=[(36,68),(73,104),(109,144),(151,177),(204,246),(250,289),(295,323),(324,335)]

cmap=plt.get_cmap('gist_rainbow').resampled(10)
colors=[[c for c in cmap(k)] for k in range(9)]
# 3D figures
for state in ['state1','state2','state3']:
    proj.chimera.close()
    data=proj['iREDbond'][f'.+AvOb.+{state}'][0]
    data.select.molsys.chimera()
    index0=np.ones(data.label.shape,dtype=bool)
    i=np.argmax(data.label==276)
    for hlx0,color in zip(hlx,colors):
        index=np.logical_and(data.label>=hlx0[0],data.label<=hlx0[1])
        index0[index]=False
        data.select.chimera(index=index,x=data.CCnorm[-2][i][index],color=color)
    data.select.chimera(index=index0,x=data.CCnorm[-2][i][index0],color=colors[-1])
        
    proj.chimera.command_line(['~ribbon','~show','show :1-337','color #1 tan',
                               '~show ~@N,C,CA,HN,O','turn x -90','view','turn y 135',
                               'show :276','color :276|:275@C,O,CA&~:276@C,O black','lighting full'])
    
    proj.chimera.savefig(os.path.join(fig_dir,f'apoCC_{state}_rho4.png'),
                         options='transparentBackground True',overwrite=True)

