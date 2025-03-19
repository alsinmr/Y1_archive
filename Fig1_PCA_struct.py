# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import os
from ClusterAna import pca
from matplotlib.colors import ListedColormap

proj=pca.project
fig_dir='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/Figure1'

#%% First, plot all clusters as histograms

def CMAP(color):   #Make a colormap that fades to white and transparent from a base color
    x=.8
    x0=1
    color=np.array(color)[:3]
    colors=[[*(np.ones(3)*x0-(x0-color)*(m/255)**0.5).tolist(),m/256*x+(1-x)] for m in range(256)]
    colors[0][-1]=0
    return ListedColormap(colors)


fig,ax=plt.subplots(2,4,sharex=True,sharey=True,figsize=[10.3,5])

PCs=[(1,2),(1,2),(3,4),(1,2)]

cmap=plt.get_cmap('tab10')

bins=np.linspace(-130,130,66)

for pca0,(pc00,pc10),ax0 in zip(pca.pcas,PCs,ax.T):
    for k in range(pca0.Cluster.n_clusters):
        i=pca0.Cluster.state==k
        for pc0,pc1,a in zip([0,pc00],[1,pc10],ax0):
            a.hist2d(pca0.PCamp[pc0][i],pca0.PCamp[pc1][i],bins=bins,cmap=CMAP(cmap(k)))
    for k in range(pca0.Cluster.n_clusters):
        for pc0,pc1,a in zip([0,pc00],[1,pc10],ax0):
            a.scatter(pca0.PCamp[pc0][pca0.Cluster.cluster_index[k]],
                           pca0.PCamp[pc1][pca0.Cluster.cluster_index[k]],
                           color='black',s=10,marker='^')
            # a.text(pca0.PCamp[pc0][pca0.Cluster.cluster_index[k]],
            #                pca0.PCamp[pc1][pca0.Cluster.cluster_index[k]],
            #                f'{pca0.Cluster.populations[k]*100:.0f}%',
            #                horizontalalignment='left' if k%2 else 'right')
            a.set_xlim([-130,130])
            a.set_ylim([-130,130])
            a.set_xlabel(f' PC {pc0}')
            a.set_ylabel(f' PC {pc1}')
            a.set_aspect('equal')
fig.tight_layout()
fig.savefig(os.path.join(fig_dir,'PCA.png'))

#%% Evaluate RMS among clusters for each trajectory
RMS=[]
for pca0 in pca.pcas:
    pos=np.array([pca0.pos[k] for k in pca0.Cluster.cluster_index])
    RMS.append(np.sqrt(((pos-pos.mean(0))**2).sum(-1).mean(0)))
    

RMS=np.array(RMS)
RMS=RMS.reshape([4,3,RMS.shape[1]//3]).mean(1)
RMS=np.log(RMS)
RMSsc=1.3*(RMS-RMS.min())/(RMS.max()-RMS.min())
RMSsc[RMSsc>1]=1
    

#%% New chimera plots
cmap=plt.get_cmap('tab20')

proj.chimera.current=0
proj.chimera.close()
for q,(pca0,chain,title,rms) in enumerate(zip(pca.pcas,['A','B','R','D'],['apo','NPY','UR-MK299','Gi'],RMSsc)):
    # if title!='apo':continue
    for k in range(pca0.Cluster.n_clusters):
        # if k:continue
        
        first=k==0 and title=='apo'
        
        pca0.traj[pca0.Cluster.cluster_index[k]]
        
        filename=os.path.join(fig_dir,'PDB.pdb')
        
        pca0.uni.atoms.write(filename)
        
        
        color=np.array([1,0,0])
        color_str=','.join([f'{100*c:.0f}' for c in color])
        # color0=np.array([210,180,140])/256
        color0=np.array([240,240,240])/256
        
        if first:
            proj.chimera.command_line([f'open {filename}',f'color #1 {color_str}','~show',
                                       'turn x -90','turn y 145','turn z 5',
                                       'alias cylinders  cartoon style  protein modeh tube rad 2 sides 24',
                                       'cylinders','lighting simple','graphics silhouettes true',
                                       'set bgColor white','lighting shadows true','view',
                                       'show #1:276&~@H*','color #1:276 black target ba'])
        else:
            proj.chimera.command_line([f'open {filename}',f'color #2 {color_str}',
                                       f'align #2/{chain}:40-290@CA toAtoms #1/A:40-290@CA',
                                       'style stick','~show','cylinders','~ribbon #1','ribbon #2',
                                       f'show #2/{chain}:276&~@H*','color #2:276 black target ba'])
            
        if title=='Gi':
            proj.chimera.command_line(['transparency #2/A,B,C 75 target r',
                                       '~ribbon #2/E','show #2/E@N,C,CA','style stick'])
        elif title=='NPY':
            proj.chimera.command_line(['~ribbon #2/A','show #2/A@N,C,CA','style stick'])
            
            
        # rms=rms.reshape([3,len(rms)//3]).mean(0) 
        cmds=[]
        for resi,rms0 in zip(pca0.atoms.residues,rms):
            c0=','.join([str(int(c*100)) for c in color*rms0+color0*(1-rms0)])
            cmds.append(f'color #{1 if first else 2}/{chain}:{resi.resid} {c0}')
        proj.chimera.command_line(cmds)
        proj.chimera.command_line(f'color #{1 if first else 2}:276 black target ba')
        
        
        filename=os.path.join(fig_dir,f'n{title}_cluster{k}.png')
        # proj.chimera.savefig(os.path.join(fig_dir,filename),options='transparentBackground True',
        #                       overwrite=True)
        proj.chimera.command_line('close #2')


    
    