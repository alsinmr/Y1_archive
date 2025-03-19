# -*- coding: utf-8 -*-


import os
from ClusterAna import pca,ECC
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pyDR


proj=pca.project

# Figure directory
fig_dir='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/Figure6'


#%% Plot time vs. rotameric state

lengths=[26359,26307,15492]
cl=np.concatenate([[0],np.cumsum(lengths)])

fig,ax=plt.subplots(1,3,sharex=True,figsize=[9,3])
i=np.argmax(ECC.resids==276)

color=[plt.get_cmap('turbo').resampled(5)(k) for k in range(1,5)]
cmapTrp=ListedColormap(color)

for k,a in enumerate(ax):
    t=ECC.t[:lengths[k]]/1e3
    a.plot(t,ECC.state[i][cl[k]:cl[k+1]],color='grey',linestyle=':')
    for m in range(4):
        q=ECC.state[i][cl[k]:cl[k+1]]==m
        a.scatter(t[q],ECC.state[i][cl[k]:cl[k+1]][q],color=cmapTrp(m),s=5)
        print(f'state {m}: {q.sum()/lengths[k]*100:.1f} %')
fig.tight_layout()


#%% Figure: Apo dynamics

proj=pyDR.Project('Projects/apo')

# apo amplitude
hlx=[(36,68),(74,102),(110,144),(154,175),(209,241),(257,289),(295,323)]
index=np.zeros(proj[-1].label.shape,dtype=bool)
for hlx0 in hlx:
    i=np.logical_and(proj[-1].label>=hlx0[0],proj[-1].label<=hlx0[1])
    index[i]=True
    


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
    
# apo CC
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