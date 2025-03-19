# -*- coding: utf-8 -*-


import os
from ClusterAna import ECC,pca
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pyDR
from Entropy import ECCnpy

proj=pca.project

# Figure directory
fig_dir='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/FigureParts/NPY'

color=[plt.get_cmap('turbo').resampled(5)(k) for k in range(1,5)]
cmapTrp=ListedColormap(color)



step=1
i=np.argmax(ECCnpy.resids==276)
pc_state_ind=[]

start=.1

from matplotlib.colors import ListedColormap

pc_amp=[(-68,-60,-86),(-20,-14,6)]

vis_frames=[]


bins=np.linspace(-130,130,50)
fig,ax=plt.subplots(2,2,figsize=[6,6],sharex=True,sharey=True)
for (n,ax0,(pc0,pc1,pc2)) in zip(range(4),ax,pc_amp):
    for n0,n1,a in zip([0,1],[1,2],ax0):
        color0=cmapTrp(n)
        colors=[[*[(1-c0)*k/256*(1-start)+c0 for c0 in color0[:-1]],1] for k in np.arange(256)]
        colors[-1][-1]=0
        cmap=ListedColormap(colors[::-1])
        q=ECCnpy.state[i]==n
        print(f'State {n}: {q.sum()/len(q)*100:.1f} %')
        # a.scatter(pca.pcas[0].PCamp[n0][q][::step],pca.pcas[0].PCamp[n1][q][::step],color=cmap(n+1),s=.1)
        a.hist2d(pca.pcas[1].PCamp[n0][q],pca.pcas[1].PCamp[n1][q],cmap=cmap,bins=bins)
        

        

        
        m=np.argmin((pca.pcas[1].PCamp[0][q]-pc0)**2+(pca.pcas[1].PCamp[1][q]-pc1)**2+\
            (pca.pcas[1].PCamp[2][q]-pc2)**2)
        m=np.argwhere(q)[:,0][m]
        if n0==0:vis_frames.append(m)
        x,y=pca.pcas[1].PCamp[[n0,n1]][:,m]
        a.scatter(x,y,color=cmap(255),marker=['o','^','*','s'][n],s=40,edgecolors='black')
        
        # q=ECC.state[i]==n
        # print(f'State {n}: {q.sum()/len(q)*100:.1f} %')
        # # a.scatter(pca.pcas[0].PCamp[n0][q][::step],pca.pcas[0].PCamp[n1][q][::step],color=cmap(n+1),s=.1)
        # a.hist2d(pca.pcas[0].PCamp[n0][q],pca.pcas[0].PCamp[n1][q],cmap=cmap,bins=bins)
        
        a.set_xlim([-130,130])
        a.set_ylim([-130,130])
        a.set_xlabel(f'PC {n0}')
        a.set_ylabel(f'PC {n1}')
            
    
fig.tight_layout()
fig.savefig(os.path.join(fig_dir,'Trp_state_vs_PC.png'))


fig,ax=plt.subplots(2,2,figsize=[6,6],sharex=True,sharey=True)
for (n,ax0,(pc0,pc1,pc2)) in zip(range(4),ax,pc_amp):
    for n0,n1,a in zip([0,1],[1,2],ax0):
        color0=cmapTrp(n)
        colors=[[*[(1-c0)*k/256*(1-start)+c0 for c0 in color0[:-1]],(k-255)/256] for k in np.arange(256)]
        colors[-1][-1]=0
        cmap=ListedColormap(colors[::-1])

        i=np.argmax(ECC.resids==276)
        q=ECC.state[i]==n
        print(f'State {n}: {q.sum()/len(q)*100:.1f} %')
        # a.scatter(pca.pcas[0].PCamp[n0][q][::step],pca.pcas[0].PCamp[n1][q][::step],color=cmap(n+1),s=.1)
        a.hist2d(pca.pcas[0].PCamp[n0][q],pca.pcas[0].PCamp[n1][q],cmap=cmap,bins=bins)

        i=np.argmax(ECCnpy.resids==276)
        q=ECCnpy.state[i]==n
        print(f'State {n}: {q.sum()/len(q)*100:.1f} %')
        # a.scatter(pca.pcas[0].PCamp[n0][q][::step],pca.pcas[0].PCamp[n1][q][::step],color=cmap(n+1),s=.1)
        
        color0=[0,0,0,1]
        colors=[[*[(1-c0)*k/256*(1-start)+c0 for c0 in color0[:-1]],(k-255)/256] for k in np.arange(256)]
        colors[-1][-1]=0
        cmap=ListedColormap(colors[::-1])
        a.hist2d(pca.pcas[1].PCamp[n0][q],pca.pcas[1].PCamp[n1][q],cmap=cmap,bins=bins)
        

        

        
        a.set_xlim([-130,130])
        a.set_ylim([-130,130])
        a.set_xlabel(f'PC {n0}')
        a.set_ylabel(f'PC {n1}')
            
    
fig.tight_layout()
fig.savefig(os.path.join(fig_dir,'Trp_state_vs_PC.png'))

#%% W 276 Time dependence

i=np.argmax(ECCnpy.resids==276)

# Plot W276 state vs. time
fig,ax=plt.subplots(1,3)
l=ECCnpy.traj.lengths
for k,a in enumerate(ax):
    state=ECCnpy.state[i,l[:k].sum():l[:k+1].sum()]

    
    a.plot(np.arange(state.size),state,color='black',linewidth=.5,linestyle=':')
    for m in range(2):
        q=state==m
        a.scatter(np.arange(state.size)[q],np.ones(q.sum())*m,s=15,color=cmapTrp(m))
        print(f'{q.sum()/q.__len__()*100:.1f}')
    a.set_xlabel('t / ns')
    a.set_yticklabels([])
    a.set_ylabel('state')
    a.set_ylim([-.1,3.1])
    a.set_xlim([-1300,27600])
fig.set_size_inches([8.7,2.9])
fig.tight_layout()
fig.savefig(os.path.join(fig_dir,'Trp276_state.png'))



#%% 3D aplitude
hlx=[(36,68),(74,102),(110,144),(154,175),(209,241),(257,289),(300,323)]

proj=pyDR.Project('Projects/npy_bound/')


cmap=plt.get_cmap('gist_earth').resampled(10)
colors=[[c for c in cmap(k)] for k in range(8)]
maxR=0
# 3D figures

from copy import copy
select=copy(pca.pcas[1].select)
select.select_bond('15N')

sub=proj['.+AvOb']['MD']

proj.chimera.current=0
for k,data in enumerate(sub):
    name=data.source.additional_info.split('_')[1]
    
    data.select=select
    data.select.molsys.make_pdb(pca.pcas[1].Cluster.cluster_index[k])
    for rho in range(4,5):
        proj.chimera.close()
        data.select.molsys.chimera()
        for hlx0,color in zip(hlx,colors):
            label=np.array([int(lbl.split('_')[1]) if lbl.split('_')[0]=='B' else 0 for lbl in data.label])
            index=np.logical_and(label>=hlx0[0],label<=hlx0[1])
            data.select.chimera(index=index,x=data.R[index][:,rho]*10,color=color)
            maxR=max(data.R[index][:,-2].max(),maxR)
            
        proj.chimera.command_line(['~ribbon','~show','show :1-337','color #1 tan',
                                   '~show ~@N,C,CA,HN,O','turn x -90','view','turn y 135',
                                   'show :276','color :276&~@N,C,CA,HN,O black',
                                   'lighting full'])
        
        
        proj.chimera.savefig(os.path.join(fig_dir,f'BB_{name}_rho{rho}.png'),
                             options='transparentBackground True',overwrite=True)
        
        
        
#%% 3D CC
hlx=[(36,68),(73,104),(109,144),(151,177),(204,246),(250,289),(295,323),(324,340)]

sub=proj['.+AvOb']['iREDbond']

cmap=plt.get_cmap('gist_rainbow').resampled(10)
colors=[[c for c in cmap(k)] for k in range(9)]
# 3D figures
from copy import copy
select=copy(pca.pcas[1].select)
select.select_bond('15N')


for k,data in enumerate(sub):
    name=data.source.additional_info.split('_')[1]
    data.select=select
    data.select.molsys.make_pdb(pca.pcas[1].Cluster.cluster_index[k])
    for rho in range(4,5):
        proj.chimera.close()
        data.select.molsys.chimera()
        index0=np.ones(data.label.shape,dtype=bool)
        i=np.argmax(data.label=='B_276')
        for hlx0,color in zip(hlx,colors):
            label=np.array([int(lbl.split('_')[1]) if lbl.split('_')[0]=='B' else 0 for lbl in data.label])
            index=np.logical_and(label>=hlx0[0],label<=hlx0[1])
            index0[index]=False
            data.select.chimera(index=index,x=data.CCnorm[rho][i][index],color=color)
        data.select.chimera(index=index0,x=data.CCnorm[rho][i][index0],color=colors[-1])
            
        proj.chimera.command_line(['~ribbon','~show','show :1-337','color #1 tan',
                                   '~show ~@N,C,CA,HN,O','turn x -90','view','turn y 135',
                                   'show :276','color :276|:275@C,O,CA&~:276@C,O black','lighting full'])
        
        proj.chimera.savefig(os.path.join(fig_dir,f'CC_{name}_rho{rho}.png'),
                             options='transparentBackground True',overwrite=True)
        
        
