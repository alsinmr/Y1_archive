# -*- coding: utf-8 -*-

from ClusterAna import pca
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import os

fig_dir='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/FigureParts/DeactivatedStates'

x=.8
x0=1

def CMAP(color):   #Make a colormap that fades to white and transparent from a base color
    color=np.array(color)[:3]
    colors=[[*(np.ones(3)*x0-(x0-color)*(m/255)**0.5).tolist(),m/256*x+(1-x)] for m in range(256)]
    colors[0][-1]=0
    return ListedColormap(colors)


fig,ax=plt.subplots()

pca.pcas[0].Hist.plot(cmap=CMAP([1,0,0]),ax=ax)
pca.pcas[2].Hist.plot(cmap=CMAP([0,0,1]),ax=ax)
ax.set_xlim([-130,130])
ax.set_ylim([-130,130])
ax.set_aspect('equal')

fig.tight_layout()
fig.savefig(os.path.join(fig_dir,'PCA.png'))

PCapo=[86,-86]
PCanta=[19.5,58.5]

ax.scatter(*PCapo,s=100,color=[1,0,0],edgecolor='black')
ax.scatter(*PCanta,s=100,color=[0,0,1],edgecolor='black')

proj=pca.project

proj.chimera.current=0
proj.chimera.close()

i=np.argmin((pca.pcas[0].PCamp[0]-PCapo[0])**2-(pca.pcas[0].PCamp[1]-PCapo[1])**2)
pca.pcas[0].traj[i]
filename=os.path.join(fig_dir,'apo.pdb')
pca.pcas[0].uni.atoms.write(filename)


proj.chimera.command_line([f'open {filename}','color #1 firebrick'])


i=np.argmin((pca.pcas[2].PCamp[0]-PCanta[0])**2-(pca.pcas[2].PCamp[1]-PCanta[1])**2)
pca.pcas[2].traj[i]
filename=os.path.join(fig_dir,'anta.pdb')
pca.pcas[2].uni.atoms.write(filename)

proj.chimera.command_line([f'open {filename}','color #2 royal blue',
                           'align #2@CA toAtoms #1@CA','turn x -90',
                           'turn y 140','view','~show','show :276&~@H*'])

proj.chimera.savefig(os.path.join(fig_dir,'view1.png'),
                     options='transparentBackground True',overwrite=True)


