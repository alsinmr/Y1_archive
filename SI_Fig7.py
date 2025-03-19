# -*- coding: utf-8 -*-


import pyDR
import matplotlib.pyplot as plt
import numpy as np
import os
from PCA import All as pca


fig_dir='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/FigureParts/DynamicsCompare'


#%% 3D plots: Amplitude analysis

hlx=[(36,68),(74,102),(110,144),(154,175),(209,241),(257,289),(300,323)]

proj=pyDR.Project('Projects/full_length_ana')


cmap=plt.get_cmap('gist_earth').resampled(10)
colors=[[c for c in cmap(k)] for k in range(8)]
maxR=0
# 3D figures


sub=proj['.+AvOb']['MD']

proj.chimera.current=0
for data,chain in zip(sub,[None,'B',None,'D']):
    name=data.source.additional_info.split('_')[1]
    for rho in range(4,5):
        proj.chimera.close()
        data.select.molsys.chimera()
        for hlx0,color in zip(hlx,colors):
            if chain is None:
                index=np.logical_and(data.label>=hlx0[0],data.label<=hlx0[1])
            else:
                label=np.array([int(lbl.split('_')[1]) if lbl.split('_')[0]==chain else 0 for lbl in data.label])
                index=np.logical_and(label>=hlx0[0],label<=hlx0[1])
            data.select.chimera(index=index,x=data.R[index][:,rho]*10,color=color)
            maxR=max(data.R[index][:,-2].max(),maxR)
            
        proj.chimera.command_line(['~ribbon','~show','show :1-337','color #1 tan',
                                   '~show ~@N,C,CA,HN,O','turn x -90','view','turn y 135',
                                   'show :276','color :276&~@N,C,CA,HN,O black',
                                   'lighting full'])
        
        if chain=='D':
            proj.chimera.command_line(['view /D,E','zoom .95'])
            
        
        
        proj.chimera.savefig(os.path.join(fig_dir,f'BB_{name}_rho{rho}.png'),
                             options='transparentBackground True',overwrite=True)
        
        
hlx=[(36,68),(73,104),(109,144),(151,177),(204,246),(250,289),(300,323),(324,335)]

cmap=plt.get_cmap('gist_earth').resampled(10)
colors=[[c for c in cmap(k)] for k in range(8)]
maxR=0
# 3D figures


sub=proj['.+AvOb']['MD']

proj.chimera.current=0
for data,chain in zip(sub,[None,'B',None,'D']):
    #if chain!='D':continue
    name=data.source.additional_info.split('_')[1]
    for rho in range(4,5):
        proj.chimera.close()
        data.select.molsys.chimera()
        index0=np.ones(data.label.shape,dtype=bool)
        for hlx0,color in zip(hlx,colors):
            if chain is None:
                index=np.logical_and(data.label>=hlx0[0],data.label<=hlx0[1])
            else:
                label=np.array([int(lbl.split('_')[1]) if lbl.split('_')[0]==chain else 0 for lbl in data.label])
                index=np.logical_and(label>=hlx0[0],label<=hlx0[1])
            index0[index]=False
            data.select.chimera(index=index,x=data.R[index][:,rho]*3,color=color)
            maxR=max(data.R[index][:,-2].max(),maxR)
        data.select.chimera(index=index0,x=data.R[index0][:,rho]*3,color=colors[-1])
            
        proj.chimera.command_line(['~ribbon','~show','show :1-337','color #1 tan',
                                   '~show ~@N,C,CA,HN,O','turn x -90','view','turn y 135',
                                   'show :276','color :276&~@N,C,CA,HN,O black',
                                   'lighting full'])
        
        if chain=='D':
            proj.chimera.command_line(['view /D,E','zoom .95'])
        
        
        proj.chimera.savefig(os.path.join(fig_dir,f'Loop_{name}_rho{rho}.png'),
                             options='transparentBackground True',overwrite=True)
        
#%% Cross-correlation analysis

hlx=[(36,68),(73,104),(109,144),(151,177),(204,246),(250,289),(295,323),(324,335)]

sub=proj['.+AvOb']['iREDbond']

cmap=plt.get_cmap('gist_rainbow').resampled(10)
colors=[[c for c in cmap(k)] for k in range(9)]
# 3D figures
for data,chain in zip(sub,[None,'B',None,'D']):
    if chain!='D':continue    
    name=data.source.additional_info.split('_')[1]
    for rho in range(4,5):
        proj.chimera.close()
        data.select.molsys.chimera()
        index0=np.ones(data.label.shape,dtype=bool)
        i=np.argmax(data.label==(276 if chain is None else f'{chain}_276')) 
        for hlx0,color in zip(hlx,colors):
            if chain is None:
                index=np.logical_and(data.label>=hlx0[0],data.label<=hlx0[1])
            else:
                label=np.array([int(lbl.split('_')[1]) if lbl.split('_')[0]==chain else 0 for lbl in data.label])
                index=np.logical_and(label>=hlx0[0],label<=hlx0[1])
            index0[index]=False
            data.select.chimera(index=index,x=data.CCnorm[rho][i][index],color=color)
        data.select.chimera(index=index0,x=data.CCnorm[rho][i][index0],color=colors[-1])
            
        proj.chimera.command_line(['~ribbon','~show','show :1-337','color #1 tan',
                                   '~show ~@N,C,CA,HN,O','turn x -90','view','turn y 135',
                                   'show :276','color :276|:275@C,O,CA&~:276@C,O black','lighting full'])
        
        if chain=='D':
            proj.chimera.command_line(['view /D,E','zoom .95'])
        
        proj.chimera.savefig(os.path.join(fig_dir,f'CC_{name}_rho{rho}.png'),
                             options='transparentBackground True',overwrite=True)
        


