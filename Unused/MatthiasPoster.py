# -*- coding: utf-8 -*-




import os
from PCA import All as pca
import pyDR
from pyDR.Entropy import EntropyCC
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


#%% Now, perform PCA and clustering of apo trajectory

proj=pyDR.Project()

md_dir='/Volumes/My Book/Y1/apo'

topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k+1}.xtc') for k in range(3)]

commands=['turn x -90','turn y 180','view','sel :276']

select=pyDR.MolSelect(topo,trajs,project=proj)
select.select_bond('15N')
select.traj.step=1


ECC=EntropyCC(select)

index=np.argmax(ECC.resids==276)

N=ECC.index[3:,index].sum()

chi=[ECC.chi[k][ECC.index[k+3,:index].sum()] for k in range(N-1,-1,-1)]
from sklearn import cluster
from sklearn.cluster import MiniBatchKMeans
cluster=MiniBatchKMeans(n_clusters=4,random_state=3)
state=cluster.fit(np.array(chi).T).labels_
ECC.state[index]=state

# Backbone clustering    
pca.pcas[0].Cluster.index=[0,1,2,3]
pca.pcas[0].Cluster.n_clusters=6
pca.pcas[0].Cluster.cluster_kwargs['random_state']=40


ax=plt.subplots()[1]

ax.cla()

# cm=plt.get_cmap('Greys')
# colors=np.array([cm(k) for k in range(255)])
# colors[:,-1]=np.linspace(0,1,255)**.25
# cm=ListedColormap(colors)

# pca.pcas[0].Hist.plot(ax=ax,cmap=cm)

# cm=plt.get_cmap('Oranges')
# colors=np.array([cm(k) for k in range(255)])
# colors[:,-1]=np.linspace(0,1,255)**.25
# cm=ListedColormap(colors)

# pca.pcas[1].Hist.plot(ax=ax,cmap=cm)

ax.scatter(pca.pcas[1].PCamp[0],pca.pcas[1].PCamp[1],s=.01,color=plt.get_cmap('Oranges').resampled(4)(2))

# cm=plt.get_cmap('Greens')
# colors=np.array([cm(k) for k in range(255)])
# colors[:,-1]=np.linspace(0,1,255)**.25
# cm=ListedColormap(colors)

# pca.pcas[2].Hist.plot(ax=ax,cmap=cm)

ax.scatter(pca.pcas[2].PCamp[0],pca.pcas[2].PCamp[1],s=.01,color=plt.get_cmap('Greens').resampled(4)(2))
ax.set_xlabel(r'PC 0 / $\AA$')
ax.set_ylabel(r'PC 1 / $\AA$')

i=np.argmax(ECC.resids==276)
cmap=plt.get_cmap('Blues').resampled(4)
colors=[[52,144,211,255],[115,192,239,255],[175,226,246,255],[0,0,0,255]]
colors=[[c for c in cmap(k)] for k in [3,2,1]]
colors.append([0,0,0,1])
# colors=[np.array(c)/256 for c in colors]
cm=ListedColormap(colors)

for k,state in enumerate([4,0,1,5]):
    index=state==ECC.state[i]
    ax.scatter(pca.pcas[0].PCamp[0][index],pca.pcas[0].PCamp[1][index],
               s=.01*(k+1),color=cm(k))
    
    
#%% Four plots
fig,ax=plt.subplots(2,4)
# ax=ax.flatten()

cmap=plt.get_cmap('Blues').resampled(4)
colors=[[52,144,211,255],[115,192,239,255],[175,226,246,255],[0,0,0,255]]
colors=[[c/256 for c in color] for color in colors]
# colors=[[c for c in cmap(k)] for k in [3,2,1]]
# colors.append([0,0,0,1])
# colors=[np.array(c)/256 for c in colors]
cm=ListedColormap(colors)

for q,ax0 in enumerate(ax):
    for k,(a,state) in enumerate(zip(ax0,[0,1,2,3])):
        a.scatter(pca.pcas[1].PCamp[1],pca.pcas[1].PCamp[2*q],s=.01,color=plt.get_cmap('Oranges').resampled(4)(2))
        
        a.scatter(pca.pcas[2].PCamp[1],pca.pcas[2].PCamp[2*q],s=.01,color=plt.get_cmap('Greens').resampled(4)(2))
        a.set_xlabel(r'PC 1 / $\AA$')
        a.set_ylabel(fr'PC {2*q} / $\AA$')
        
        i=np.argmax(ECC.resids==276)
            
    
        index=state==ECC.state[i]
        a.scatter(pca.pcas[0].PCamp[1][index],pca.pcas[0].PCamp[2*q][index],
                   s=.01*(k+1),color=cm(k))
        a.set_xlim([-120,120])
        a.set_ylim([-120,120])
fig.set_size_inches([18.5,8.4])
fig.tight_layout()

    
#%% 3D scatterplot

fig=plt.figure()
ax=fig.add_subplot(1,1,1,projection='3d')

ax.cla()
step=10
ax.scatter(pca.pcas[1].PCamp[0][::step],pca.pcas[1].PCamp[1][::step],pca.pcas[1].PCamp[2][::step],
           s=.1,color=plt.get_cmap('Oranges').resampled(4)(2))

ax.scatter(pca.pcas[2].PCamp[0][::step],pca.pcas[2].PCamp[1][::step],pca.pcas[2].PCamp[2][::step],
           s=.1,color=plt.get_cmap('Greens').resampled(4)(2))

i=np.argmax(ECC.resids==276)
cmap=plt.get_cmap('Blues').resampled(4)
colors=[[c for c in cmap(k)] for k in [3,2,1]]
colors.append([0,0,0,1])
# colors=[np.array(c)/256 for c in colors]
cm=ListedColormap(colors)

for k,state in enumerate([4,0,1,5]):
    index=state==ECC.state[i]
    ax.scatter3D(pca.pcas[0].PCamp[0][index][::step],pca.pcas[0].PCamp[1][index][::step],
                 pca.pcas[0].PCamp[2][index][::step],
               s=.1*(k+1),color=cm(k))


    
#%% 3D plots
# Figure directory
proj=pyDR.Project('Projects/apo')

fig_dir='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Report240828/figure_parts/apo'

hlx=[(36,68),(74,102),(110,144),(154,175),(209,240),(257,289),(295,323)]

cmap=plt.get_cmap('gist_earth').resampled(10)
colors=[[c for c in cmap(k)] for k in range(8)]
maxR=0
# 3D figures
for state in ['state1','state2','state3']:
    proj.chimera.close()
    data=proj['MD'][f'.+AvOb.+{state}'][0]
    data.select.molsys.chimera()
    for hlx0,color in zip(hlx,colors):
        index=np.logical_and(data.label>=hlx0[0],data.label<=hlx0[1])
        data.select.chimera(index=index,x=data.R[index][:,-2]*10,color=color)
        maxR=max(data.R[index][:,-2].max(),maxR)
        
    proj.chimera.command_line(['~ribbon','~show','show :1-337','color #1 tan',
                               '~show ~@N,C,CA,HN,O','turn x -90','view','turn y 135',
                               'show :276','color :276&~@N,C,CA,HN,O black',
                               'lighting full'])
    
    proj.chimera.savefig(os.path.join(fig_dir,f'apo_{state}_rho4_new.png'),
                         options='transparentBackground True',overwrite=True)