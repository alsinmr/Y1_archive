# -*- coding: utf-8 -*-

import pyDR
import os
from pyDR.PCA.PCAclean import PCA

from pyDR.Entropy import EntropyCC

import matplotlib.pyplot as plt
import numpy as np




TRP_res=[106,148,163,276,288]

proj=pyDR.Project('Projects/PCA_3state/',create=True)

md_dir='/Volumes/My Book/Y1/apo'

topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k+1}.xtc') for k in range(3)]

# Get length of each trajectory
l_traj=[]
for traj in trajs:
    ms=pyDR.MolSys(topo,traj)
    l_traj.append(len(ms.traj))


commands=['turn x -90','turn y 180','view','sel :276']

folder='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/FigureParts' #Location for  figures

select=pyDR.MolSelect(topo,trajs,project=proj)
select.select_bond('15N')
select.traj.step=1

ECC=EntropyCC(select)
ECC.load()


pca=PCA(select)
pca.select_atoms('name N C CA')
nc=6
pca.Cluster.n_clusters=nc
pca.Cluster.index=[0,1,2]
fig=pca.Cluster.plot(percent=True)[0].figure
fig.set_size_inches([8,3.6])
fig.tight_layout()

ECC.PCA=pca

transitions=[(1,2),(0,1,3),(3,4),(4,5)]

for transition in transitions:
    ECC.CCchimera(indexCC='pca',states=transition,scaling=2)
ECC.CCchimera(indexCC='pca',scaling=2)
