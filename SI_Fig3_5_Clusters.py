# -*- coding: utf-8 -*-

from ClusterAna import pca,lim,n_clusters,RMS,RMS_W276
import matplotlib.pyplot as plt
import numpy as np




def ClusterPlot(cmap='tab10'):
    fig,ax=plt.subplots(4,5,sharex=True,sharey=True,figsize=[12.1,9.5])
    ax=ax.T
    maxbin=np.abs(pca.PCamp).max()
    for k,a in enumerate(ax[0]):
        pca.Hist.plot(n0=k,n1=k+1,ax=a,maxbin=maxbin)
    
    if isinstance(cmap,str):cmap=plt.get_cmap(cmap)
    
    titles=['Apo','NPY bound','UR-MK299 bound','NPY/Gi bound']
    
    
    for k,(ax0,pca0,title) in enumerate(zip(ax[1:],pca.pcas,titles)):
        for m,a in enumerate(ax0):
            pca.Hist.plot(n0=m,n1=m+1,cmap='Greys',ax=a,maxbin=maxbin)
            for state in np.arange(pca0.Cluster.n_clusters):
                i=pca0.Cluster.state==state
                
                a.scatter(pca0.PCamp[m][i],pca0.PCamp[m+1][i],s=.01,color=cmap(state))
                
            a.set_xlim(-lim,lim)
            a.set_ylim(-lim,lim)
            a.set_xlabel(fr'PC {m} / $\AA$')
            a.set_ylabel(fr'PC {m+1} / $\AA$')
        
        ax0[0].set_title(title)
    
    # fig.set_size_inches()
    fig.tight_layout()
    
    return fig




def RMSplot(cmap='turbo_r',cmap_traj='tab10',ax=None):
    if ax is None:ax=plt.subplots()[1]
    if isinstance(cmap,str):cmap=plt.get_cmap(cmap)
    if isinstance(cmap_traj,str):cmap_traj=plt.get_cmap(cmap_traj)
    
    hdl=ax.imshow(RMS,cmap=cmap)
    plt.colorbar(hdl)
    
    ax.set_xticks(np.arange(0,21))
    ax.set_yticks(np.arange(0,21))
    
    lw=4
    
    for k in range(4):
        start=np.sum(n_clusters[:k])-.5
        stop=np.sum(n_clusters[:k+1])-.5
        
        ax.plot([start,start],[start,stop],color=cmap_traj(k),linewidth=lw)
        ax.plot([start,stop],[start,start],color=cmap_traj(k),linewidth=lw)
        ax.plot([stop,stop],[start,stop],color=cmap_traj(k),linewidth=lw)
        ax.plot([start,stop],[stop,stop],color=cmap_traj(k),linewidth=lw)
        
    return ax


# def RMS_W276plot(cmap='tab10'):
#     fig,ax=plt.subplots(4,1,figsize=[7.2,9.8],sharex=True,sharey=True)
#     labels=[]
#     if isinstance(cmap,str):cmap=plt.get_cmap(cmap)
#     for k,title in enumerate(['apo','NPY','UR-MK299','Gi']):
#         for q in range(n_clusters[k]):
#             labels.append(f'{title} {q}')
#     for q,(rms,a) in enumerate(zip(RMS_W276,ax)):
#         a.set_title(f'W276 state {q}')
#         for k in range(4):
#             i0=np.sum(n_clusters[:k]).astype(int)
#             i1=np.sum(n_clusters[:k+1]).astype(int)
#             a.bar(np.arange(i0,i1),rms[i0:i1],color=cmap(k))
#         a.set_xticks(np.arange(np.sum(n_clusters)))
#         a.set_xticklabels(labels,rotation=90)
#         a.set_ylabel(r'RMS / $\AA$')
#     fig.tight_layout()
    
#     return fig
    
# def RMS_W276plot(cmap='tab10'):
#     fig,ax=plt.subplots(4,1,figsize=[7.2,9.8],sharey=True)
    
#     if isinstance(cmap,str):cmap=plt.get_cmap(cmap)
    
#     for k,(title,a) in enumerate(zip(['apo','NPY','UR-MK299','Gi'],ax)):
#         n=n_clusters[k]
#         i0=np.sum(n_clusters[:k]).astype(int)
#         i1=np.sum(n_clusters[:k+1]).astype(int)
#         label=[]
#         for q,rms in enumerate(RMS_W276):
#             a.bar(np.arange(n)+q*n,rms[i0:i1],color=cmap(q))
#             label.extend([f'clstr {m}' for m in range(n)])
#         a.set_xticks(np.arange(n*4))
#         a.set_xticklabels(label,rotation=90)
#         a.set_title(title)
#     ax[0].legend([f'W276 {k+1}' for k in range(4)])
        
#     fig.tight_layout()
#     return fig

def RMS_W276plot(cmap='tab10'):
    fig,ax=plt.subplots(4,1,figsize=[7.2,9.8],sharey=True)
    
    if isinstance(cmap,str):cmap=plt.get_cmap(cmap)
    
    for k,(title,a) in enumerate(zip(['apo','NPY','UR-MK299','Gi'],ax)):
        n=n_clusters[k]
        start=np.sum(n_clusters[:k]).astype(int)
        label=[]
        for q in range(n):
            a.bar(np.arange(4)+q*4,RMS_W276[:,start+q],color=cmap(q))
            label.extend([f'W276 {m}' for m in range(4)])
        a.legend([f'clstr {m}' for m in range(n)])
                  
        a.set_xticks(np.arange(len(label)))
        a.set_xticklabels(label,rotation=90)
        a.set_title(title)
        a.set_ylabel(r'RMS / $\AA$')
    ax[0].legend([f'W276 {k+1}' for k in range(4)])
        
    fig.tight_layout()
    return fig
    
            