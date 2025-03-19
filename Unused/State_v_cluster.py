# -*- coding: utf-8 -*-

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

x=.8
x0=1

def CMAP(color):   #Make a colormap that fades to white and transparent from a base color
    color=np.array(color)[:3]
    colors=[[*(np.ones(3)*x0-(x0-color)*(m/255)**0.5).tolist(),m/256*x+(1-x)] for m in range(256)]
    colors[0][-1]=0
    return ListedColormap(colors)
    
#%% State and cluster

bins=np.linspace(-130,130,66)

i=np.argmax(ECC.resids==276)
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