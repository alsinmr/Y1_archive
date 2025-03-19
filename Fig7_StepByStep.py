# -*- coding: utf-8 -*-

# Make a plot going through the activation states of Y1

import os
from ClusterAna import pca
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pyDR
from Entropy import ECC,ECCnpy

ECC=[ECC,ECCnpy,None,None]


# Step 1 : Y1 apo deactivated, trp state 1
# Step 2 : Y1 apo near-activated, trp state 2
# Step 3 : Y1 NPY-bound, near activated, trp state 2
# Step 4 : Y1 NPY-bound, activated, trp state 1
# Step 5 : Y1 Gi-bound, activated


# Subparts : Structures, W276 orientation
# Histograms
# Detector analysis
# CC to W276 (??)

fig_dir='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/Figure7/'

proj=pca.project

# List of PC positions, trajectory, and trp state


PCs=[(84,-84),(28,-20),(-20,-12),(-44,16),(-71,-35)]   #Positions in PCA to select frames
Wstates=[0,1,1,0,0]                   #States of W276 to select from
flips=[False,True,True,False,False]   #Does W276 flip for this transition?
trajs=[0,0,1,1,3]                     #Trajectory to select from

cmap0=plt.get_cmap('jet').resampled(7)  #Build a color map with desired colors
colors=[cmap0(k) for k in [0,2,3,5,6]]
colors[-1]=plt.get_cmap('jet').resampled(20)(18)
cmap0=ListedColormap(colors)

bins=np.linspace(-130,130,66) #Binning for histograms





x=.8
x0=1

def CMAP(color):   #Make a colormap that fades to white and transparent from a base color
    color=np.array(color)[:3]
    colors=[[*(np.ones(3)*x0-(x0-color)*(m/255)**0.5).tolist(),m/256*x+(1-x)] for m in range(256)]
    colors[0][-1]=0
    return ListedColormap(colors)



#%% Plot all frames where flipping occurs
ax=plt.subplots()[1]
pca0=pca.pcas[0]
ECC0=ECC[0]
i=np.argmax(ECC0.resids==276)     #Find W276
q=ECC0.state[i]==0 # Frames with correct state of W276
PCamp0,PCamp1=pca0.PCamp[0][q],pca0.PCamp[1][q] #PC amplitudes
cmap=CMAP(cmap0(0))    #Get the color map
ax.hist2d(PCamp0,PCamp1,cmap=cmap,bins=bins)    #Histogram

q=ECC0.state[i]==1 # Frames with correct state of W276
PCamp0,PCamp1=pca0.PCamp[0][q],pca0.PCamp[1][q] #PC amplitudes
cmap=CMAP(cmap0(1))    #Get the color map
ax.hist2d(PCamp0,PCamp1,cmap=cmap,bins=bins)    #Histogram


frames=np.argwhere(np.logical_and(ECC0.state[i][:-1]==0,ECC0.state[i][1:]==1))[:,0]
ax.scatter(pca0.PCamp[0][frames+1],pca0.PCamp[1][frames+1])
    

ax=plt.subplots()[1]
pca0=pca.pcas[1]
ECC0=ECC[1]
i=np.argmax(ECC0.resids==276)     #Find W276
q=ECC0.state[i]==0 # Frames with correct state of W276
PCamp0,PCamp1=pca0.PCamp[0][q],pca0.PCamp[1][q] #PC amplitudes
cmap=CMAP(cmap0(0))    #Get the color map
ax.hist2d(PCamp0,PCamp1,cmap=cmap,bins=bins)    #Histogram

q=ECC0.state[i]==1 # Frames with correct state of W276
PCamp0,PCamp1=pca0.PCamp[0][q],pca0.PCamp[1][q] #PC amplitudes
cmap=CMAP(cmap0(1))    #Get the color map
ax.hist2d(PCamp0,PCamp1,cmap=cmap,bins=bins)    #Histogram


frames=np.argwhere(np.logical_and(ECC0.state[i][:-1]==1,ECC0.state[i][1:]==0))[:,0]
for frame in frames:
    if frame==16735:continue
    pc0,pc1=pca0.PCamp[0][frame],pca0.PCamp[1][frame]
    ax.scatter(pc0,pc1,s=20,color=cmap0(0),edgecolor='black')
    


#%% Make histograms, save structures
s=100
edgecolor='black'
lw=2

p=[]
for p0 in enumerate(zip(PCs,Wstates,trajs,flips)):
    p.append(p0)
    
p=[p[0],p[1],p[3],p[2],p[4]]
# p=[p[0],p[1]]



fig,ax=plt.subplots(1,4,sharex=True,sharey=True,figsize=[16.3,4.3])
for k,((pc0,pc1),Wstate,traj,flip) in p:
    pca0=pca.pcas[traj]    #Pick the correct PCA
    ECC0=ECC[traj]         #Pick the correct W267 state 
    cmap=CMAP(cmap0(k))    #Get the color map
    
    if ECC0 is not None:
        i=np.argmax(ECC0.resids==276)     #Find W276
        q=ECC0.state[i]==Wstate # Frames with correct state of W276
        PCamp0,PCamp1=pca0.PCamp[0][q],pca0.PCamp[1][q] #PC amplitudes
        print(f'Step {k+1}: W276 population for state {Wstate} is {q.sum()/len(q)*100:.1f} %')
    else:
        PCamp0,PCamp1=pca0.PCamp[0],pca0.PCamp[1]  #PC amplitudes (no W276 state selection)
    
    
    if flip:
        if k==1:
            frames=np.argwhere(np.logical_and(ECC0.state[i][:-1]==0,ECC0.state[i][1:]==1))[:,0]
            m=np.argmin((pca0.PCamp[0][frames+1]-pc0)**2+(pca0.PCamp[1][frames+1]-pc1)**2)
            frame=frames[m]+1
        else:
            frames=np.argwhere(np.logical_and(ECC0.state[i][:-1]==1,ECC0.state[i][1:]==0))[:,0]
            frames=frames[frames!=16735]
            m=np.argmin((pca0.PCamp[0][frames]-pc0)**2+(pca0.PCamp[1][frames]-pc1)**2)
            # m=5
            frame=frames[m]
        
    else:
        m=np.argmin((PCamp0-pc0)**2+(PCamp1-pc1)**2)
        frame=np.argwhere(q)[:,0][m] if ECC0 is not None else m
    pc0,pc1=pca0.PCamp[0][frame],pca0.PCamp[1][frame]
    
    # if k<len(ax):
    #     ax[k].hist2d(PCamp0,PCamp1,cmap=cmap,bins=bins)    #Histogram
    #     ax[k].scatter(pc0,pc1,s=s,color=cmap0(k),edgecolor='black')
    #     if flip:
    #         ax[k].scatter(pca0.PCamp[0][frame+1],pca0.PCamp[1][frame+1],s=s,color='grey',edgecolor='black')
    

    if k==0:
        ax[k].hist2d(PCamp0,PCamp1,cmap=cmap,bins=bins)
        ax[k].scatter(pc0,pc1,s=s,color=cmap0(k),edgecolor=edgecolor,linewidth=lw)
    elif k==1:

        for a in ax[k-1:k+1]:
            a.hist2d(PCamp0,PCamp1,cmap=cmap,bins=bins)
            a.scatter(pc0,pc1,s=s,color=cmap0(k),edgecolor=edgecolor,linewidth=lw) 
        ax[k-1].scatter(pca0.PCamp[0][frame-1],pca0.PCamp[1][frame-1],s=s,color='black',marker='x',edgecolor=edgecolor,linewidth=lw)
    elif k==2:
        for a in ax[k-1:k+1]:
            a.hist2d(PCamp0,PCamp1,cmap=cmap,bins=bins)
            a.scatter(pc0,pc1,s=s,color=cmap0(k),edgecolor=edgecolor,linewidth=lw)
        
        ax[k].scatter(pca0.PCamp[0][frame+1],pca0.PCamp[1][frame+1],s=s,color='black',edgecolor=edgecolor,marker='x',linewidth=lw)
    else:
        for a in ax[k-1:k+1]:
            a.hist2d(PCamp0,PCamp1,cmap=cmap,bins=bins)
            a.scatter(pc0,pc1,s=s,color=cmap0(k),edgecolor=edgecolor,linewidth=lw)
            
    # Write out selected pdb
    pca0.traj[frame]
    filename=os.path.join(fig_dir,f'Y1_step{k}.pdb')
    pca0.uni.atoms.write(filename)
    # Also write out next frame if flip occurs
    if flip:
        pca0.traj[frame+(-1 if k==1 else 1)]
        filename=os.path.join(fig_dir,f'Y1_step{k}_flip.pdb')
        pca0.uni.atoms.write(filename)
        
    
    
ax[0].set_ylabel('PC 1')
for a in ax:a.set_xlabel('PC 0')

fig.tight_layout()
fig.savefig(os.path.join(fig_dir,'PCA_hist.png'))


# Plot in chimeraX
proj.chimera.current=0
proj.chimera.close()
proj.chimera.command_line(['~show'])
for k,segid in enumerate(['A','A','B','B','D']):
    if k:proj.chimera.command_line(['~show','~ribbon'])
    filename=os.path.join(fig_dir,f'Y1_step{k}.pdb')
    color=','.join([f'{100*c:.0f}' for c in cmap0(k)[:3]])
    
    
    
    
    if k:
        proj.chimera.command_line([f'open {filename}',f'color #2 {color}','~show #2',
                                   f'show #2/{segid}:276&~@H*',f'color #20/{segid}:276&~@H* black target ba',
                                   'ribbon #2'])
        proj.chimera.command_line(f'align #2/{segid}:28-337@CA toAtoms #1/A:28-337@CA')
    else:
        proj.chimera.command_line([f'open {filename}',f'color #1 {color}','~show #1',
                                   f'show #1/{segid}:276&~@H*',f'color #10/{segid}:276&~@H* black target ba',
                                   'ribbon #1','lighting full'])
        proj.chimera.command_line(['turn x -90','turn y 125','turn z 0','view','move z -45'])
        
    if k in [1,2]:
        filename=os.path.join(fig_dir,f'Y1_step{k}_flip.pdb')
        proj.chimera.command_line([f'open {filename}','color #3 grey','~show #3',
                                   f'show #3/{segid}:276&~@H*',f'color #3/{segid}:276&~@H* grey target ba',
                                   'ribbon #3'])
        proj.chimera.command_line(f'align #3/{segid}:28-337@CA toAtoms #1/A:28-337@CA')
        

                               
    if k==4:
        # proj.chimera.command_line([f'transparency #{k+1}&~/D,E 50 target r','view',
        #                            'style stick','zoom 1.7','move y -25'])
        proj.chimera.command_line(['transparency #2&~/D,E 50 target r',
                                   'style stick'])
        
    proj.chimera.savefig(os.path.join(fig_dir,f'Y1_BB_step{k}.png'),
                         options='transparentBackground True',overwrite=True)
    
    proj.chimera.command_line(['turn y 70','zoom 3','move x -5'])
    
    proj.chimera.savefig(os.path.join(fig_dir,f'Y1_W276_step{k}.png'),
                         options='transparentBackground True',overwrite=True)
    
    proj.chimera.command_line(['move x 5','zoom 0.3333','turn y -70','close #2','close #3'])
    
    
# State change/ ring flip time dependence
fig,ax=plt.subplots(2,3,sharex=True)

pca0=pca.pcas[0]

l=pca0.traj.lengths
for k,a in enumerate(ax[0]):
    state=pca0.Cluster.state[np.sum(l[:k]):np.sum(l[:k+1])]
    a.scatter(np.arange(len(state)),state,s=.1,color='black')
    a.set_ylim([-.1,5.1])
    
i=np.argmax(ECC[0].resids==276)
for k,a in enumerate(ax[1]):
    state=ECC[0].state[i][np.sum(l[:k]):np.sum(l[:k+1])]
    a.scatter(np.arange(len(state)),state,s=.1,color='black')
    a.set_ylim([-.1,3.1])
    
    
fig,ax=plt.subplots(2,3,sharex=True)

pca0=pca.pcas[1]

l=pca0.traj.lengths
for k,a in enumerate(ax[0]):
    state=pca0.Cluster.state[np.sum(l[:k]):np.sum(l[:k+1])]
    a.scatter(np.arange(len(state)),state,s=.1,color='black')
    a.set_ylim([-.1,5.1])
    
i=np.argmax(ECC[1].resids==276)
for k,a in enumerate(ax[1]):
    state=ECC[1].state[i][np.sum(l[:k]):np.sum(l[:k+1])]
    a.scatter(np.arange(len(state)),state,s=.1,color='black')
    a.set_ylim([-.1,3.1])
    
    
    

    

