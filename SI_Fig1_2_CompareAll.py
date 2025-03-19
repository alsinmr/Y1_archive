# -*- coding: utf-8 -*-

import os
from PCA import All as pca
import pyDR
import numpy as np
import matplotlib.pyplot as plt


# Figure directory
fig_dir='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Report240828/figure_parts/all'

#%% Figure: Comparing trajs
proj=pca.project

fig=pca.hist_by_traj(nmax=6)
fig.set_size_inches([10.3,13])
fig.tight_layout()
fig.savefig(os.path.join(fig_dir,'pca_all_0to6.png'))

fig=pca.hist_by_traj(nmax=1)
fig.set_size_inches([11.5,2.7])
fig.tight_layout()

# Example structures
# PCs=[(52.5,-20,22.5,28),(-43.5,-55,-81,-.5),(20.6,57.6,-21.8,-24.4),(-68,-37,9.2,-20.8)]

PCs=[(86,-86),(-42.5,14.5),(20,56),(-71,-35.5)]

cmap=plt.get_cmap('tab10')
markers=['*','^','p','P']

proj.chimera.current=0
proj.chimera.close()

for k,(a,pca0,PC) in enumerate(zip(fig.axes,pca.pcas,PCs)):
    error=np.sum([(PC0-PCamp)**2 for PC0,PCamp in zip(PC,pca0.PCamp)],axis=0)
    i=np.argmin(error)
    a.scatter(pca0.PCamp[0][i],pca0.PCamp[1][i],s=100,edgecolor='black',color=cmap(k),marker=markers[k])
    
    pca0.traj[i]
    pca0.uni.atoms.write(os.path.join(fig_dir,'temp.pdb'))
    
    color=','.join([f'{int(c*100)}' for c in cmap(k)[:-1]])
    proj.chimera.command_line([f'open {os.path.join(fig_dir,"temp.pdb")}',
                               f'color #{k+1} {color}'])
    
proj.chimera.command_line(['style stick','~show','ribbon','~ribbon #4&~/D,E',
                           'align #2/B:28-339@CA toAtoms #1:28-339@CA',
                           'align #3/R:28-339@CA toAtoms #1:28-339@CA',
                           'align #4/D:28-339@CA toAtoms #1:28-339@CA',
                           'show :276&~@H*','~show #4&~/D,E',
                           'turn x -90','turn y 135','turn z -5','view'])

proj.chimera.savefig(os.path.join(fig_dir,'Y1_view0.png'),options='transparentBackground True',overwrite=True)
proj.chimera.command_line(['turn y 180'])
proj.chimera.savefig(os.path.join(fig_dir,'Y1_view1.png'),options='transparentBackground True',overwrite=True)
proj.chimera.command_line(['turn y -88','zoom 5','move y 3','transparency 50 target r'])
proj.chimera.savefig(os.path.join(fig_dir,'Y1_view2.png'),options='transparentBackground True',overwrite=True)
proj.chimera.command_line(['transparency 0 target r','move y -3','zoom 0.2','turn y 88'])

for k in range(2,4):proj.chimera.command_line(['align #{k}:28-339@CA toAtoms #1@CA:28-339'])
proj.chimera.command_line(['align #4/D,E@CA toAtoms #1@CA'])
      
fig.savefig(os.path.join(fig_dir,'pca_all.png'))
