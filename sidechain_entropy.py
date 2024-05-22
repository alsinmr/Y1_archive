#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 16:42:13 2024

@author: albertsmith
"""

import pyDR
import os
import matplotlib.pyplot as plt
import numpy as np
from copy import copy
from pyDR.misc.tools import AA

from pyDR.Entropy import EntropyCC,CombineEntropy

TRP_res=[106,148,163,276,288]

proj=pyDR.Project()

md_dir='/Volumes/My Book/Y1/apo'

topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k+1}.xtc') for k in range(3)]

commands=['turn x -90','turn y 180','view','sel :276']

folder='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/sidechain_ana' #Location for  figures

select=pyDR.MolSelect(topo,trajs,project=proj)
select.select_bond('15N')
select.traj.step=1

if os.path.exists('EntropyCC/CC_apo.data'):
    ECC=EntropyCC('EntropyCC/CC_apo.data')
    ECC.project=proj
else:
    ECC=EntropyCC(select)
    ECC.save('EntropyCC/CC_apo.data')




fig,ax=plt.subplots()
ECC.plotCC(ax=ax)
fig.savefig(os.path.join(folder,'apo_totalCC.png'))




proj.chimera.close()
ECC.chimera(norm=False)
proj.chimera.command_line(commands)
proj.chimera.savefig(os.path.join(folder,'Y1_entropy.png'),options='transparentBackground true')

proj.chimera.close()
ECC.chimera(norm=True)
proj.chimera.command_line(commands)
proj.chimera.savefig(os.path.join(folder,'Y1_entropy_norm.png'),options='transparentBackground true')

for trp in TRP_res:
    proj.chimera.close()
    ECC.CCchimera(indexCC=np.argmax(ECC.resids==trp))
    proj.chimera.command_line(commands)
    proj.chimera.savefig(os.path.join(folder,f'Y1_{trp}CC.png'),options='transparentBackground true',overwrite=True)

#%% Plot time dependence of largest correlations
i0=np.argmax(ECC.resids==276)
i=np.argsort(ECC.CC[i0])[-4:-1]

fig,ax=plt.subplots(2,2)
ax=ax.flatten()
ax[0].plot(ECC.state[i0])
ax[0].text(4000,.3,'W276')
ax[0].set_ylabel('State')

lines=[725,1025,1400,19175,
       26359,26359+432,26359+8700,26359+19705,26359+21000,
       52666]

yl=ax[0].get_ylim()
ax[0].set_ylim(yl)
for line in lines:
    ax[0].plot(np.ones(2)*line,yl,color='grey',linestyle=':')

for ii,ax0 in zip(i,ax[1:]):
    ax0.plot(ECC.state[ii])
    ax0.text(4000,.3,f'{AA(ECC.resi[ii].resname).symbol}{ECC.resids[ii]}, CC={ECC.CC[i0,ii]:.2f}')
    if ax0.is_last_row():
        ax0.set_xlabel('Frame')
    if ax0.is_first_col():
        ax0.set_ylabel('State')
fig.savefig(os.path.join(folder,'Y1_276_time_depend.png'))

#%% Separate according to W276 state
from States import chunks

if os.path.exists('EntropyCC/CC_state1.data'):
    ECCc=[EntropyCC(f'EntropyCC/CC_state{k+1}.data') for k in range(3)]
    for ECC0 in ECCc:ECC0.project=proj
else:
    ECCs=[[] for _ in range(3)]
    select=[pyDR.MolSelect(topo,traj,project=proj).select_bond('15N') for traj in trajs]
    
    
    length=6600
    for run,start,state,count in zip(*chunks(length).values()):
        print([run,start,state,count])
        ECCs[state-1].append(EntropyCC(select[run]))
        ECCs[state-1][-1].select.traj.t0=start
        ECCs[state-1][-1].select.traj.tf=start+length
        _=ECCs[state-1][-1].vt
        
    ECCc=[]
    for k,ECC0 in enumerate(ECCs):
        ECCc.append(CombineEntropy(*ECC0))
        ECCc[-1].save(f'EntropyCC/CC_state{k+1}.data')
    
for k,ECC0 in enumerate(ECCc):
    proj.chimera.close()
    ECC0.CCchimera(indexCC=np.argmax(ECC.resids==276))
    proj.chimera.command_line(commands)
    proj.chimera.savefig(os.path.join(folder,f'Y1_CC_state{k+1}.png'),options='transparentBackground true',overwrite=True)
    
for k,ECC0 in enumerate(ECCc):
    fig,ax=plt.subplots()
    ECC0.plotCC(ax=ax)
    fig.savefig(os.path.join(folder,f'state{k+1}_CC.png'))
    
    CC=copy(ECC0.CC)
    CC[np.isnan(CC)]=0
    
    res_max=ECC0.resids[np.argmax(CC.sum(0))]
    
    proj.chimera.close()
    ECC0.CCchimera(indexCC=np.argmax(ECC.resids==res_max))
    proj.chimera.command_line(commands)
    proj.chimera.savefig(os.path.join(folder,f'Y1_CC_state{k+1}_{res_max}.png'),options='transparentBackground true',overwrite=True)
    
#%% NPY bound
if os.path.exists('EntropyCC/CC_NPY.data'):
    ECCnpy=EntropyCC('EntropyCC/CC_NPY.data')
    ECCnpy.project=proj
else:
    md_dir='/Volumes/My Book/Y1/NPY'
    
    topo=os.path.join(md_dir,'Y1R_NPY.pdb')
    trajs=[os.path.join(md_dir,f'run{k+1}.proc.xtc') for k in range(3)]
    
    select=pyDR.MolSelect(topo,trajs,project=proj)
    select.select_bond('15N')
    select.traj.step=1
    
    ECCnpy=EntropyCC(select)
    ECCnpy.save('EntropyCC/CC_NPY.data')


fig,ax=plt.subplots()
ECCnpy.plotCC(ax=ax)
fig.savefig(os.path.join(folder,'NPY_CC.png'))

CC=copy(ECCnpy.CC)
CC[np.isnan(CC)]=0

res_max=ECCnpy.resids[np.argmax(CC.sum(0))]

proj.chimera.close()
ECCnpy.CCchimera(indexCC=np.argmax(ECCnpy.resids==res_max))
proj.chimera.command_line(commands)
proj.chimera.savefig(os.path.join(folder,f'Y1_CC_NPY_{res_max}.png'),options='transparentBackground true',overwrite=True)

#%% NPY/Gprotein bound
if os.path.exists('EntropyCC/CC_Gi.data'):
    ECCgi=EntropyCC('EntropyCC/CC_Gi.data')
    ECCgi.project=proj
else:
    md_dir='/Volumes/My Book/Y1/NPY_Gi'
    
    topo=os.path.join(md_dir,'Y1R_NPY_Gi.pdb')
    trajs=[os.path.join(md_dir,f'run{k+1}.proc.xtc') for k in range(3)]
    
    select=pyDR.MolSelect(topo,trajs,project=proj)
    select.select_bond('15N',segids=['D','E'])
    select.traj.step=5
    
    ECCgi=EntropyCC(select)
    # ECCgi.save('EntropyCC/CC_Gi.data')

fig,ax=plt.subplots()
ECCgi.plotCC(ax=ax)
fig.savefig(os.path.join(folder,'NPY_gi_CC.png'))

CC=copy(ECCgi.CC)
CC[np.isnan(CC)]=0

res_max=ECCgi.resids[np.argmax(CC.sum(0))]

proj.chimera.close()
ECCgi.CCchimera(indexCC=np.argmax(CC.sum(0)))
proj.chimera.command_line(commands)
proj.chimera.savefig(os.path.join(folder,f'Y1_CC_gi_{res_max}.png'),options='transparentBackground true',overwrite=True)

# Plot time dependence of R33 in NPY
i=328 #(R33)
fig,ax=plt.subplots()
ax.plot(ECCgi.state[i])
ax.set_xlabel('Frame')
ax.set_ylabel('State')
ax.set_title('NPY R33 in Gi bound NPY')
fig.savefig(os.path.join(folder,'Y1_Gi_NPY33_t_depend.png'))

i1=np.argmax(ECCgi.resids==254)
fig,ax=plt.subplots(2,1)
ax[0].plot(ECCgi.state[i])
ax[0].set_xlabel('Frame')
ax[0].set_ylabel('State')
ax[1].plot(ECCgi.state[i1])
ax[1].set_xlabel('Frame')
ax[1].set_ylabel('State')
ax[0].set_title('R33_v_R254')

fig.savefig(os.path.join(folder,'R33_v_R254_Gi_bound.png'))