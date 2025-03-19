# -*- coding: utf-8 -*-

from ClusterAna import pca
import pyDR
import numpy as np
import matplotlib.pyplot as plt



l=10000

state=np.zeros([4,3,l],dtype=np.uint8)

# Collect the states
for k,pca0 in enumerate(pca.pcas):
    for q in range(3):
        i0=np.sum(pca0.traj.lengths[:q+1])-l
        state[k,q]=pca0.Cluster.state[i0:i0+l]
        
# Calculate correlation
Ct=np.ones([4,3,l-1])

for k in range(4):
    for q in range(3):

        
        for m in range(1,l-1):
            p=np.unique(state[k,q][m:],return_counts=True)[1]/(l-m)
            S0=-(p*np.log(p)).sum()
            p=np.unique(state[k,q][:-m],return_counts=True)[1]/(l-m)
            S0+=-(p*np.log(p)).sum()
            
            state0=state[k,q,m:]*6+state[k,q,:-m]
            p=np.unique(state0,return_counts=True)[1]/(l-m)
            S=-(p*np.log(p)).sum()
            Ct[k,q,m]=2*(S0-S)/S0
            
        print([k,q])
 
    
N=np.logical_not(np.isnan(Ct)).sum(1)

Ct_avg=np.zeros([4,l-1])
for k in range(4):
    for q in range(3):
        i=np.logical_not(np.isnan(Ct[k,q]))
        Ct_avg[k,i]+=Ct[k,q,i]/N[k,i]
        
Ct_avg=Ct_avg[:,:l//2]




t=np.arange(l//2)

ax=plt.subplots()[1]
ax.plot(t,Ct_avg.T)



data=pyDR.Data(Ct_avg,sens=pyDR.Sens.MD(t=t))


data.detect.r_auto(5)
fit=data.fit()
opt=fit.opt2dist(rhoz_cleanup=True)
po=opt.plot(style='bar')
po.show_tc()