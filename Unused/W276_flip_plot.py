# -*- coding: utf-8 -*-

from ClusterAna import pca
from Entropy import ECC
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

apo=pca.pcas[0]

#%% Plot exit from particular rotameric state
color=[plt.get_cmap('turbo').resampled(5)(k) for k in range(1,5)]
cmapTrp=ListedColormap(color)
cmap=plt.get_cmap('tab20')

x=.8
x0=1

def CMAP(color):   #Make a colormap that fades to white and transparent from a base color
    color=np.array(color)[:3]
    colors=[[*(np.ones(3)*x0-(x0-color)*(m/255)**0.5).tolist(),m/256*x+(1-x)] for m in range(256)]
    colors[0][-1]=0
    return ListedColormap(colors)


bins=np.linspace(-130,130,66)
step=1
nsteps=100

fig,ax=plt.subplots(2,2,sharex=True,sharey=True,figsize=[8,8])
ax=ax.flatten()

stay=50

i=np.argmax(ECC.resids==276)
for state,a in enumerate(ax):
    rot_index=np.argwhere(np.logical_and(ECC.state[i][:-1]==state,ECC.state[i][1:]!=state))[:,0]
    q=np.logical_not(np.array([ri in apo.traj.lengths-1 for ri in rot_index])) #Exclude changes between trajectories
    rot_index=rot_index[q]
    # check if we hop back to original state
    q=np.logical_not(np.array([state in ECC.state[i][ri+1:ri+stay+1] for ri in rot_index])) #Exclude short traverses
    rot_index=rot_index[q]
    q=np.array([np.all(state==ECC.state[i][ri-stay:ri]) for ri in rot_index]) #Exclude short traverses
    rot_index=rot_index[q]
    
    index=ECC.state[i]==state
    a.hist2d(apo.PCamp[0][index],apo.PCamp[1][index],bins=bins,cmap=CMAP(cmapTrp(state)))
    
    for k,ri in enumerate(rot_index):
        a.scatter(apo.PCamp[0][ri],apo.PCamp[1][ri],color='black',marker='x',s=5)
        
        
        l=np.cumsum([0,*apo.traj.lengths.tolist()])
        ti=np.argmax(ri<l[1:])
        
        a.text(apo.PCamp[0][ri],apo.PCamp[1][ri],f'{ti},{ri-l[ti]}')
        
        start=max(ri-nsteps*step,l[ti])
        stop=min(ri+nsteps*step,l[ti+1])
        
        a.plot(apo.PCamp[0][start:stop:step],apo.PCamp[1][start:stop:step],
                   color=cmap(k%20))
    
    
#%% State and cluster
fig,ax=plt.subplots(2,2,sharex=True,sharey=True,figsize=[8,8])
ax=ax.flatten()
cmap=plt.get_cmap('tab10')
for state,a in enumerate(ax):
    index=ECC.state[i]==state
    index[:26358]=False
    index[26358+26300:]=False
    for cluster in range(6):
        q=apo.Cluster.state[index]==cluster
        frac=q.sum()/index.sum()
        a.hist2d(apo.PCamp[0][index][q],apo.PCamp[1][index][q],bins=bins,cmap=CMAP(cmap(cluster)))
        # a.scatter(apo.PCamp[0][index][q],apo.PCamp[1][index][q],color=cmap(cluster),s=.1)
        if np.any(q):
            pc0,pc1=apo.PCamp[0][index][q].mean(),apo.PCamp[1][index][q].mean()
            a.scatter(pc0,pc1,marker='^',color=cmap(cluster),s=50)
            a.text(pc0,pc1,f'{frac*100:.1f}%')
        
        # print(f'W276 state {state+1}, cluster {cluster}:{frac*100:.1f}%')