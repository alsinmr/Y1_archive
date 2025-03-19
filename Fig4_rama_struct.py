# -*- coding: utf-8 -*-

import numpy as np
from Entropy import ECC,ECCnpy
from PCA import All as pca
from PCA import step
import matplotlib.pyplot as plt
import os
import pyDR
from matplotlib.colors import ListedColormap

proj=ECC.project

color=[plt.get_cmap('turbo').resampled(5)(k) for k in [1,2,3,4,0]]
cmap=ListedColormap(color)


#%% Apo
index=np.argmax(ECC.resids==276)
N=ECC.index[3:,index].sum()
chi=[ECC.chi[k][ECC.index[k+3,:index].sum()] for k in range(N-1,-1,-1)]

fig,ax=plt.subplots(figsize=[5,4])


for state,k in zip(range(4),[0,1,2,3]):
    i=(ECC.state[index]==state).flatten()
    ax.scatter(chi[0][i]+120,chi[1][i],s=.1,c=cmap(k))
ax.set_xlim([0,360])
ax.set_ylim([0,360])
ax.set_xticks(np.linspace(0,360,7))
ax.set_yticks(np.linspace(0,360,7))
ax.set_aspect(1)
ax.set_xlabel(r'$\chi_1$ / $^\circ$')
ax.set_ylabel(r'$\chi_2$ / $^\circ$')


fig.tight_layout()

fig.savefig('/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/Figure4/apo_rama.png',transparent=True)


#%% NPY-bound



index=np.argmax(ECCnpy.resids==276)
N=ECCnpy.index[3:,index].sum()
chi=[ECCnpy.chi[k][ECCnpy.index[k+3,:index].sum()] for k in range(N-1,-1,-1)]

fig,ax=plt.subplots(figsize=[5,4])


for state,k in zip(range(4),[0,1,2,3]):
    i=(ECCnpy.state[index]==state).flatten()
    ax.scatter(chi[0][i]+120,chi[1][i],s=.1,c=cmap(k))
ax.set_xlim([0,360])
ax.set_ylim([0,360])
ax.set_xticks(np.linspace(0,360,7))
ax.set_yticks(np.linspace(0,360,7))
ax.set_aspect(1)
ax.set_xlabel(r'$\chi_1$ / $^\circ$')
ax.set_ylabel(r'$\chi_2$ / $^\circ$')


fig.tight_layout()
fig.savefig('/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/Figure4/npy_rama.png',transparent=True)

# #%% NPY/Gi ana
# from pyDR.Entropy import EntropyCC

# md_dir_as='/Volumes/My Book/Y1/NPY_Gi'

# topo=os.path.join(md_dir_as,'Y1R_NPY_Gi.pdb')
# trajs=[os.path.join(md_dir_as,f'run{k+1}.proc.xtc') for k in range(3)]

# select=pyDR.MolSelect(topo,trajs,project=proj)
# select.select_bond('15N',segids='D',resids=276)
# select.traj.step=step

# ECC_Gi=EntropyCC(select)

# ECC_Gi.PCA=pca.pcas[2]

# index=np.argmax(ECC_Gi.resids==276)
# N=ECC_Gi.index[3:,index].sum()
# chi=[ECC_Gi.chi[k][ECC_Gi.index[k+3,:index].sum()] for k in range(N-1,-1,-1)]


# from sklearn.cluster import MiniBatchKMeans
# cluster=MiniBatchKMeans(n_clusters=2,random_state=5)
# state=cluster.fit(np.array(chi).T).labels_
# # i=np.logical_and(state==3,chi[1]>325)
# # state[i]=0
# ECC_Gi.state[index]=state

# fig,ax=plt.subplots(figsize=[5,4])

# for state,k in zip(range(4),[0,2]):
#     i=(ECC_Gi.state[index]==state).flatten()
#     ax.scatter(chi[0][i]+120,chi[1][i],s=.1,c=cmap(k))
# ax.set_xlim([0,360])
# ax.set_ylim([0,360])
# ax.set_xticks(np.linspace(0,360,7))
# ax.set_yticks(np.linspace(0,360,7))
# ax.set_aspect(1)
# ax.set_xlabel(r'$\chi_1$ / $^\circ$')
# ax.set_ylabel(r'$\chi_2$ / $^\circ$')


# fig.tight_layout()
# fig.savefig('/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/Figure4/gi_rama.png',transparent=True)

# #%% antagonist ana
# from pyDR.Entropy import EntropyCC

# md_dir_as='/Volumes/My Book/Y1/Antagonist'

# topo=os.path.join(md_dir_as,'Y1R_UR_MK299.pdb')
# trajs=[os.path.join(md_dir_as,f'run{k+1}.xtc') for k in range(3)]

# select=pyDR.MolSelect(topo,trajs,project=proj)
# select.select_bond('15N',segids='R',resids=276)
# select.traj.step=step

# ECCanta=EntropyCC(select)

# ECCanta.PCA=pca.pcas[2]

# index=np.argmax(ECCanta.resids==276)
# N=ECCanta.index[3:,index].sum()
# chi=[ECCanta.chi[k][ECCanta.index[k+3,:index].sum()] for k in range(N-1,-1,-1)]


# from sklearn.cluster import MiniBatchKMeans
# cluster=MiniBatchKMeans(n_clusters=3,random_state=5)
# state=cluster.fit(np.array(chi).T).labels_
# state[:]=0
# state[chi[0]<220]=1
# state[chi[0]<100]=2
# # i=np.logical_and(state==3,chi[1]>325)
# # state[i]=0
# ECCanta.state[index]=state

# fig,ax=plt.subplots(figsize=[5,4])

# for state,k in zip(range(4),[4,0,1,3,4]):
#     i=(ECCanta.state[index]==state).flatten()
#     ax.scatter(np.mod(chi[0][i]+120,360),chi[1][i],s=.1,c=cmap(k))
# ax.set_xlim([0,360])
# ax.set_ylim([0,360])
# ax.set_xticks(np.linspace(0,360,7))
# ax.set_yticks(np.linspace(0,360,7))
# ax.set_aspect(1)
# ax.set_xlabel(r'$\chi_1$ / $^\circ$')
# ax.set_ylabel(r'$\chi_2$ / $^\circ$')


# fig.tight_layout()
# fig.savefig('/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/Figure4/anta_rama.png',transparent=True)


#%% Find representative frames

# Apo
index=np.argmax(ECC.resids==276)
N=ECC.index[3:,index].sum()
chi=[ECC.chi[k][ECC.index[k+3,:index].sum()] for k in range(N-1,-1,-1)]



proj=ECC.project
proj.chimera.current=0
proj.chimera.close()

apo_index=[]
for state in range(3):
    i=(ECC.state[index]==state).flatten()
    chi1,chi2=chi[0][i].mean(),chi[1][i].mean()
    e=(chi[0]-chi1)**2+(chi[1]-chi2)**2
    apo_index.append(np.argmin(e))
    
for i in apo_index:
    ECC.select.molsys.make_pdb(ti=i)
    ECC.select.molsys.chimera()
    
for k in range(2,4):
    proj.chimera.command_line(f'align #{k}:274-277@CA toAtoms #1:274-277@CA')

for state,k in zip(range(1,4),[0,1,2]):
    c=','.join([str(int(100*c)) for c in cmap(k)[:3]])
    proj.chimera.command_line(f'color #{state} {c} target ab')
    
proj.chimera.command_line(['~ribbon','~show','ribbon #1:256-289','show :276&~@H*',
                           'color light grey target r',
                           'turn x -90','turn y -70','sel :276','view','~sel',
                           'zoom 2'])

proj.chimera.savefig('/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/Figure4/apo_rotamers.png',
                     options='transparentBackground True',overwrite=True)


# NPY
index=np.argmax(ECCnpy.resids==276)
N=ECCnpy.index[3:,index].sum()
chi=[ECCnpy.chi[k][ECCnpy.index[k+3,:index].sum()] for k in range(N-1,-1,-1)]
i=(ECCnpy.state[index]==0).flatten()
chi1,chi2=chi[0][i].mean(),chi[1][i].mean()
e=(chi[0]-chi1)**2+(chi[1]-chi2)**2
npy_index=np.argmin(e)

filename=ECCnpy.select.molsys.make_pdb(ti=npy_index,replace=False)
proj.chimera.command_line(f'open "{filename}"')
proj.chimera.command_line('align #4/B:274-277@CA toAtoms #1:274-277@CA')
c=','.join([str(int(100*c)) for c in cmap(0)[:3]])
proj.chimera.command_line(f'color #4 {c} target ab')

proj.chimera.command_line(['~ribbon','~show','ribbon #4:256-289','show #4/B:276&~@H*',
                           'color light grey target r'])

proj.chimera.savefig('/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/Figure4/npy_rotamers.png',
                     options='transparentBackground True',overwrite=True)

# NPY with two rotamers
index=np.argmax(ECCnpy.resids==276)
N=ECCnpy.index[3:,index].sum()
chi=[ECCnpy.chi[k][ECCnpy.index[k+3,:index].sum()] for k in range(N-1,-1,-1)]
i=(ECCnpy.state[index]==1).flatten()
chi1,chi2=chi[0][i].mean(),chi[1][i].mean()
e=(chi[0]-chi1)**2+(chi[1]-chi2)**2
npy_index=np.argmin(e)

filename=ECCnpy.select.molsys.make_pdb(ti=npy_index,replace=False)
proj.chimera.command_line(f'open "{filename}"')
proj.chimera.command_line('align #5/B:274-277@CA toAtoms #1:274-277@CA')
c=','.join([str(int(100*c)) for c in cmap(1)[:3]])
proj.chimera.command_line(f'color #5 {c} target ab')

proj.chimera.command_line(['~ribbon','~show','ribbon #5:256-289','show #5/B:276&~@H*',
                           'color light grey target r'])

proj.chimera.savefig('/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/Figure4/npy_2rotamers.png',
                     options='transparentBackground True',overwrite=True)



# # Gi
# index=np.argmax(ECC_Gi.resids==276)
# N=ECC_Gi.index[3:,index].sum()
# chi=[ECC_Gi.chi[k][ECC_Gi.index[k+3,:index].sum()] for k in range(N-1,-1,-1)]
# i=(ECC_Gi.state[index]==0).flatten()
# chi1,chi2=chi[0][i].mean(),chi[1][i].mean()
# e=(chi[0]-chi1)**2+(chi[1]-chi2)**2
# Gi_index=np.argmin(e)

# filename=ECC_Gi.select.molsys.make_pdb(ti=Gi_index,replace=False)
# proj.chimera.command_line(f'open "{filename}"')
# proj.chimera.command_line('align #5/D:274-277@CA toAtoms #1:274-277@CA')
# c=','.join([str(int(100*c)) for c in cmap(0)[:3]])
# proj.chimera.command_line(f'color #5 {c} target ab')

# proj.chimera.command_line(['~ribbon','~show','ribbon #5/D:256-289','show #5/D:276&~@H*',
#                            'color light grey target r','style stick'])

# proj.chimera.savefig('/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/Figure4/Gi_rotamers.png',
#                      options='transparentBackground True',overwrite=True)

    

# anta
# index=np.argmax(ECCanta.resids==276)
# N=ECCanta.index[3:,index].sum()
# chi=[ECCanta.chi[k][ECCanta.index[k+3,:index].sum()] for k in range(N-1,-1,-1)]
# i=(ECCanta.state[index]==1).flatten()
# chi1,chi2=chi[0][i].mean(),chi[1][i].mean()
# e=(chi[0]-chi1)**2+(chi[1]-chi2)**2
# anta_index=np.argmin(e)

# filename=ECCanta.select.molsys.make_pdb(ti=anta_index,replace=False)
# proj.chimera.command_line(f'open "{filename}"')
# proj.chimera.command_line('align #6/R:274-277@CA toAtoms #1:274-277@CA')
# c=','.join([str(int(100*c)) for c in cmap(0)[:3]])
# proj.chimera.command_line(f'color #6 {c} target ab')

# proj.chimera.command_line(['~ribbon','~show','ribbon #6/R:256-289','show #6/R:276&~@H*',
#                            'color light grey target r','style stick'])
    
# proj.chimera.savefig('/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/Figure4/anta_rotamers.png',
#                      options='transparentBackground True',overwrite=True)
    
    


