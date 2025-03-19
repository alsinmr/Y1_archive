# -*- coding: utf-8 -*-

import os
from ClusterAna import ECC,pca
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pyDR


proj=pca.project

# Figure directory 

fig_dir='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/Figure5'

color=[plt.get_cmap('turbo').resampled(5)(k) for k in range(1,5)]
cmapTrp=ListedColormap(color)


#%% Plot histograms of each W276 state
i=np.argmax(ECC.resids==276)

fig,ax=plt.subplots(4,2,figsize=[6,10.5],sharex=True,sharey=True)
step=1
i=np.argmax(ECC.resids==276)
pc_state_ind=[]

start=.1

from matplotlib.colors import ListedColormap

pc_amp=[(86.5,-86,22),(28,-19,24),(62,1,-24),(55,1.5,8)]

vis_frames=[]
theta=60*np.pi/180

for (n,ax0,(pc0,pc1,pc2)) in zip(range(4),ax,pc_amp):
    for n0,n1,a in zip([0,1],[1,2],ax0):
        color0=cmapTrp(n)
        colors=[[*[(1-c0)*k/256*(1-start)+c0 for c0 in color0[:-1]],1] for k in np.arange(256)]
        colors[-1][-1]=0
        cmap=ListedColormap(colors[::-1])
        q=ECC.state[i]==n
        print(f'State {n}: {q.sum()/len(q)*100:.1f} %')
        # a.scatter(pca.pcas[0].PCamp[n0][q][::step],pca.pcas[0].PCamp[n1][q][::step],color=cmap(n+1),s=.1)
        a.hist2d(pca.pcas[0].PCamp[n0][q],pca.pcas[0].PCamp[n1][q],cmap=cmap,bins=50)
        
        a.set_xlim([-130,130])
        a.set_ylim([-130,130])
        a.set_xlabel(f'PC {n0}')
        a.set_ylabel(f'PC {n1}')
        
        m=np.argmin((pca.pcas[0].PCamp[0][q]-pc0)**2+(pca.pcas[0].PCamp[1][q]-pc1)**2+\
            (pca.pcas[0].PCamp[2][q]-pc2)**2)
        m=np.argwhere(q)[:,0][m]
        if n0==0:vis_frames.append(m)
        x,y=pca.pcas[0].PCamp[[n0,n1]][:,m]
        a.scatter(x,y,color=cmap(255),marker=['o','^','*','s'][n],s=40,edgecolors='black')
        if n0==0:
            pc01=pca.pcas[0].PCamp[n0][q]*np.cos(theta)-pca.pcas[0].PCamp[n1][q]*np.sin(theta)
            pcavg=(pca.pcas[0].PCamp[n0][q]*np.sin(theta)+pca.pcas[0].PCamp[n1][q]*np.cos(theta)).mean()
            # pc01-=pcavg
            bins=np.linspace(-150,200,101)
            A=np.histogram(pc01,bins=bins)[0]
            x0=(bins[:-1]+bins[1:])/2
            x=np.cos(-theta)*x0-np.sin(-theta)*(A/A.max()*50+pcavg)
            y=np.sin(-theta)*x0+np.cos(-theta)*(A/A.max()*50+pcavg)
            a.plot(x,y,color='grey',linestyle=':',linewidth=1)
            
            
    
fig.tight_layout()
fig.savefig(os.path.join(fig_dir,'Trp_state_vs_PC.png'))


#%% 3D plot of selected positions above

proj.chimera.current=0
proj.chimera.close()

for k,frame in enumerate(vis_frames):
    filename=os.path.join(fig_dir,f'Trp_state{k}.pdb')
    
    pca.pcas[0].traj[frame]
    pca.pcas[0].uni.atoms.write(filename)
    clr=[int(c*100) for c in cmapTrp(k)[:-1]]
    proj.chimera.command_line([f'open {filename}',f'color #{k+1} {clr[0]},{clr[1]},{clr[2]}',
                               'select ',f'contacts #{k+1}:276 restrict #{k+1}&~@CA,C,O,H*,N ignoreHiddenModels true color #000000'])
    if k:
        proj.chimera.command_line(f'align #{k+1}@CA toAtoms #1@CA')
        
proj.chimera.command_line(['transparency 30 target rba','~show','show :276&~@H*','~sel','color :276 grey target ba'])

# I manually reoriented the view in chimera

for k in range(4):
    proj.chimera.command_line(['~show','~ribbon',f'ribbon #{k+1}',f'show #{k+1}&~@H*',f'show #{k+1}.1&~@H*'])
    proj.chimera.savefig(os.path.join(fig_dir,f'Trp_state{k}.png'),options='transparentBackground True',overwrite=True)
    
