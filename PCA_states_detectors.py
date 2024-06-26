# -*- coding: utf-8 -*-

import pyDR
import os
from pyDR.PCA.PCAclean import PCA

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


pca=PCA(select)
pca.select_atoms('name N C CA')
nc=3
pca.Cluster.n_clusters=nc
pca.Cluster.index=[0,1,2]
fig=pca.Cluster.plot(percent=True)[0].figure
fig.set_size_inches([8,3.6])
fig.tight_layout()



#%% Plot the time-dependence of the PCA cluster
ax=plt.subplots()[1]

ax.set_ylim([-.5,2.5])
ax.set_xlabel(r't / $\mu$s')
for pos in np.cumsum(l_traj):
    ax.plot([pos/1e3,pos/1e3],ax.get_ylim(),color='grey')


titles=[f'PCA{k}' for k in range(3)]
starts=[35000,10000,55000]
length=13000

for k,start in enumerate(starts):
    color=[c for c in plt.get_cmap('tab10')(k)]
    color[-1]=.5
    ax.fill_between([start/1e3,(start+length)/1e3],[ax.get_ylim()[0],ax.get_ylim()[0]],[ax.get_ylim()[1],ax.get_ylim()[1]],
                    color=color)
    ax.text((start+length/2)/1e3,.5,f'PCA{k}',rotation=90)

ax.scatter(pca.t/1e3,pca.Cluster.state,s=5)

#%% Run iRED/direct simulation

for start in starts:
    select.traj.t0=start
    select.traj.tf=start+length
    pyDR.md2data(select)
    pyDR.md2iRED(select).iRED2data()
    
proj.detect.r_no_opt(10)
proj.fit()
proj['no_opt'].detect.r_auto(6)
proj['no_opt'].fit().opt2dist(rhoz_cleanup=True).modes2bonds()
proj.save()


