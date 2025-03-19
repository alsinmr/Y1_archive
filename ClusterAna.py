# -*- coding: utf-8 -*-


from PCA import All as pca

import numpy as np
from Entropy import ECC




#%% Cluster the apo state

indices=[[0,1,2,3],[0,1,2],[0,1,2,3,4],[1,2]]
n_clusters=[6,4,5,5]
random_state=[40,41,9,7]

pca.Cluster.index=[0,1,2,3,4]
pca.Cluster.n_clusters=np.sum(n_clusters)
state=np.zeros(len(pca.traj),dtype=int)
l=pca.lengths

for k,(pca0,index,nc,rs) in enumerate(zip(pca.pcas,indices,n_clusters,random_state)):
    pca0.Cluster.index=index
    pca0.Cluster.n_clusters=nc
    pca0.Cluster.cluster_kwargs['random_state']=rs
    # if k==1:
        # s0=pca0.Cluster.state
        # s=np.zeros(s0.shape,dtype=int)
        # s[s0==0]=1
        # s[s0==1]=3
        # s[s0==2]=0
        # s[s0==3]=2
        # pca0.Cluster._state=s
    if k==3:
        pca0.Cluster.output.labels_[pca0.Cluster.output.labels_==3]=0
        pca0.Cluster.output.labels_[pca0.Cluster.output.labels_>=3]-=1
        pca0.Cluster.output.n_clusters-=1
        n_clusters[-1]-=1
        
        
        
    state[l[:k].sum():l[:k+1].sum()]=pca0.Cluster.state+np.sum(n_clusters[:k])

pca.Cluster._state=state

lim=np.abs(pca.PCamp).max()*1.05



PCavg=np.array([pca.PCamp[:,pca.Cluster.state==state].mean(-1) for state in range(np.sum(n_clusters))])
Posavg=np.array([pca.pos[pca.Cluster.state==state].mean(0) for state in range(np.sum(n_clusters))])

RMS=np.array([np.sqrt(((S-PCavg)**2).mean(-1)) for S in PCavg])
RMSpos=np.array([np.sqrt(((S-Posavg)**2).mean(-1).mean(-1)) for S in Posavg])



i=np.argmax(ECC.resids==276)
RMS_W276=[]
for state in range(4):
    index=ECC.state[i]==state
    RMS_W276.append([np.sqrt(((S-pca.pcas[0].PCamp[:,index].T)**2).mean()) for S in PCavg])
RMS_W276=np.array(RMS_W276)


            
        
    
