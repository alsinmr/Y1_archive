# -*- coding: utf-8 -*-


from PCA import All as pca
from PCA import step
import os
import pyDR
from pyDR.Entropy import EntropyCC
import numpy as np


proj=pca.project


#%% Apo ana
md_dir='/Volumes/My Book/Y1/apo' 

topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k+1}.xtc') for k in range(3)]

commands=['turn x -90','turn y 180','view','sel :276']

select=pyDR.MolSelect(topo,trajs,project=proj)
select.select_bond('15N')
select.traj.step=step

ECC=EntropyCC(select)

ECC.PCA=pca.pcas[0]

index=np.argmax(ECC.resids==276)
N=ECC.index[3:,index].sum()
chi=[ECC.chi[k][ECC.index[k+3,:index].sum()] for k in range(N-1,-1,-1)]

from sklearn.cluster import MiniBatchKMeans
cluster=MiniBatchKMeans(n_clusters=4,random_state=3)
state=cluster.fit(np.array(chi).T).labels_
i=np.logical_and(state==3,chi[1]>325)
state[i]=0
ECC.state[index]=state

#%% NPY ana

md_dir_as='/Volumes/My Book/Y1/NPY'

topo=os.path.join(md_dir_as,'Y1R_NPY.pdb')
trajs=[os.path.join(md_dir_as,f'run{k+1}.proc.xtc') for k in range(3)]

select=pyDR.MolSelect(topo,trajs,project=proj)
select.traj.step=step

select=pyDR.MolSelect(topo,trajs,project=proj)
select.select_bond('15N')
select.traj.step=step

ECCnpy=EntropyCC(select)

ECCnpy.PCA=pca.pcas[1]

index=np.argmax(ECCnpy.resids==276)
N=ECCnpy.index[3:,index].sum()
chi=[ECCnpy.chi[k][ECCnpy.index[k+3,:index].sum()] for k in range(N-1,-1,-1)]
chi[1]=np.mod(chi[1]+180,360)


from sklearn.cluster import MiniBatchKMeans
cluster=MiniBatchKMeans(n_clusters=2,random_state=5)
state=cluster.fit(np.array(chi).T).labels_
i=np.logical_and(state==3,chi[1]>325)
state[i]=0
ECCnpy.state[index]=state


