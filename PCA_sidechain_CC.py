#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:20:52 2024

@author: albertsmith
"""


import pyDR
import os
from pyDR.PCA.PCAclean import PCA

import matplotlib.pyplot as plt
import numpy as np
from copy import copy
from pyDR.misc.tools import AA

from pyDR.Entropy import EntropyCC,CombineEntropy

TRP_res=[106,148,163,276,288]

proj=pyDR.Project('Projects/sidechain/')

md_dir='/Volumes/My Book/Y1/apo'

topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k+1}.xtc') for k in range(3)]

commands=['turn x -90','turn y 180','view','sel :276']

folder='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Report240603' #Location for  figures

select=pyDR.MolSelect(topo,trajs,project=proj)
select.select_bond('15N')
select.traj.step=1


ECC=EntropyCC(select)
ECC.load()


ECC.PCA=PCA(select)
ECC.PCA.select_atoms('name N C CA')
nc=6
ECC.PCA.Cluster.n_clusters=nc
ECC.PCA.Cluster.index=[0,1,2]
fig=ECC.PCA.Cluster.plot()[0].figure
fig.set_size_inches([8,3.6])
fig.tight_layout()
fig.savefig(os.path.join(folder,'BackbonePCA_cluster.png'))


#%% Plot the clusters
ECC.PCA.Cluster.plot()


def getW276state(i):
    from States import states
    states=np.concatenate([run for run in states.values()])
    return states[i]
    

W276state=getW276state(ECC.PCA.Cluster.cluster_index)

# Bar plot PCA state vs. Trp state
ax=plt.subplots()[1]
nt=ECC.PCA.Cluster.state.size
p=np.array([[(getW276state(ECC.PCA.Cluster.state==k)==q).sum()/nt for q in range(4)] for k in range(nc)])
hdl=[]
hatch=['','xxxx','////','oooo']
for k in range(4):
    for q in range(nc):
        color=[c for c in plt.get_cmap('tab10')(q)]
        # color[-1]=.25+k/4*.75
        h0=ax.bar(q,p[q,k:].sum()*100,color=color,hatch=hatch[k])
    hdl.append(h0)
ax.set_xlabel('PCA state')
ax.legend(hdl,('Unassigned','State 1','State 2','State 3'))
ax.set_title('W276 state vs. PCA state')
ax.set_ylabel('%')
fig=ax.figure
fig.set_size_inches([5.15,4.2])
fig.tight_layout()
fig.savefig(os.path.join(folder,'PCA_v_W276_percent.png'))



#%% Plot Cross-correlation
ax=ECC.plotCCpca()
CCpca=ECC.CCpca
a=np.argsort(CCpca)[-6:]
for a0 in a:
    ax.plot([a0,a0],CCpca[a0]+np.array([.01,.03]),color='red')
    ax.text(a0,CCpca[a0]+.05,f'{AA(ECC.resi[a0].resname).symbol}{ECC.resids[a0]}',
            rotation=90,horizontalalignment='center',verticalalignment='center')
fig=ax.figure
fig.set_size_inches([12.2,4.0])
fig.tight_layout()
fig.savefig(os.path.join(folder,'PCA_CC.png'))
    
#%% Sidechain CC vs. PCA CC
fig,ax=plt.subplots()
ax.plot(ECC.CC.sum(1)/ECC.CC.sum(1).max()*ECC.CCpca.max(),color='black')
ECC.plotCCpca(color='red',linestyle=':',ax=ax)
ax.legend(('Sum of C.C.','C.C. vs. PCA'))
ax.set_ylim(ax.get_ylim())
fig.set_size_inches([10,4])
fig.tight_layout()
fig.savefig(os.path.join(folder,'CC_PCA_vs_CCsum.png'))


fig=ECC.plotCC().figure
fig.savefig(os.path.join(folder,'sidechain_CC.png'))


#%% Plots for each state
state=   [0,1, 2,    3,   4,5]
state276=[1,2,[2,3],[1,2],1,1]

d0=6

commands=['show :276&~@H*','transparency 50 target r','turn z 135','turn x 10','view','zoom 1.5']
for s0,s2760 in zip(state,state276):
    proj.chimera.close()
    for s00 in np.atleast_1d(s2760):
        d=np.sqrt(((ECC.PCA.PCamp[ECC.PCA.Cluster.index].T-ECC.PCA.Cluster.PCavg[s0])**2).sum(-1))
        i=np.logical_and(np.logical_and(d<d0,getW276state(np.arange(nt))==s00),
                                   ECC.PCA.Cluster.state==s0)
        print(np.argwhere(i)[0,0])
        ECC.PCA.Cluster.chimera(frame=np.argwhere(i)[0,0])
    proj.chimera.command_line(commands)
    proj.chimera.savefig(os.path.join(folder,f'PCA_state{s0}_EC.png'),options='transparentBackground True',overwrite=True)
    proj.chimera.command_line('turn x 170')
    proj.chimera.savefig(os.path.join(folder,f'PCA_state{s0}_IC.png'),options='transparentBackground True',overwrite=True)
    
    proj.chimera.command_line(['turn x 80','turn y 40','view','zoom 1.1'])
    proj.chimera.savefig(os.path.join(folder,f'PCA_state{s0}.png'),options='transparentBackground True',overwrite=True)
    
#%% Time trajectory
from States import states
l=[0]
for v in states.values():l.append(v.size)
l=np.cumsum(l)
fig,ax=plt.subplots(1,3)
cmap=plt.get_cmap('turbo')
step=10
for k,a in enumerate(ax):
    t=np.arange(l[k+1]-l[k])*ECC.traj.dt
    state=getW276state(np.arange(l[k],l[k+1]))[::step]
    c=np.linspace(0,1,l[k+1]-l[k])[::step]
    ECC.PCA.Hist.plot(cmap='binary',ax=a)
    for q,m in zip(range(4),['.','x','o','D']):
        i=state==q
        a.scatter(ECC.PCA.PCamp[0,l[k]:l[k+1]:step][i],ECC.PCA.PCamp[1,l[k]:l[k+1]:step][i],
                  c=c[i],cmap=plt.get_cmap('turbo'),s=20,marker=m)
    # a.set_xlim([-150,150])
    # a.set_ylim([-150,150])
    
    
#%% Overlap with activated state
md_dir_as='/Volumes/My Book/Y1/NPY'

topo=os.path.join(md_dir_as,'Y1R_NPY.pdb')
trajs=[os.path.join(md_dir_as,f'run{k+1}.proc.xtc') for k in range(3)]

select=pyDR.MolSelect(topo,trajs,project=proj)
select.select_bond('15N')
select.traj.step=1

AS=PCA(select)
AS.select_atoms('name N C CA')

i=np.argmin(((AS.pos-AS.mean)**2).sum(-1).sum(-1)) #The most "average" frame in the bound state

ref=AS.pos[i,AS.atoms.segids=='B']

error=[]
for k,i in enumerate(ECC.PCA.Cluster.cluster_index):
    if k==1:i=616
    ECC.PCA.traj[i]
    pos=AS.align(ref,AS.atoms[AS.atoms.segids=='B'],ECC.PCA.atoms[40:])
    error.append(((ref-pos)**2).sum())

#%% Compare state overlap
commands=['align #1/B:30-300@CA toAtoms #2:30-300@CA','~show',
                           'ribbon','show :276&~@H*','~ribbon #1/A','view','lighting full']
for state in range(6):
    proj.chimera.close()
    AS.Hist.chimera(PCamp=[0])
    proj.chimera.command_line('color grey')
    
    ECC.PCA.Cluster.chimera(state=state)
    
    proj.chimera.command_line(commands)
    proj.chimera.savefig(f'/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Report240603/inactive{state}_v_activ_v2.png',
                      options='transparentBackground True',overwrite=True)
    proj.chimera.command_line(['turn y 90','turn z 90','view'])
    proj.chimera.savefig(f'/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Report240603/inactive{state}_v_activ.png',
                      options='transparentBackground True',overwrite=True)

i=616
proj.chimera.close()
AS.Hist.chimera(PCamp=[0])
proj.chimera.command_line('color grey')

ECC.PCA.Hist.chimera(PCamp=ECC.PCA.PCamp[:,i])
proj.chimera.command_line(commands)
proj.chimera.savefig('/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Report240603/inactive1_v_activ_v2.png',
                  options='transparentBackground True',overwrite=True)

proj.chimera.command_line(['turn y 90','turn z 90','view'])
proj.chimera.savefig('/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Report240603/inactive1_v_activ.png',
                  options='transparentBackground True',overwrite=True)


#%% Activation parameter
# I don't really see that it's so activated...
d26=[]
TM6=[252,267]
TM2=[75,88]
a2=ECC.PCA.atoms[np.logical_and(ECC.PCA.atoms.resids>=TM2[0],ECC.PCA.atoms.resids<=TM2[1])]
a6=ECC.PCA.atoms[np.logical_and(ECC.PCA.atoms.resids>=TM6[0],ECC.PCA.atoms.resids<=TM6[1])]
for k,i in enumerate(ECC.PCA.Cluster.cluster_index):
    if k==1:i=616
    ECC.PCA.traj[i]
    
    d26.append(np.sqrt(((a2.positions.mean(0)-a6.positions.mean(0))**2).sum()))



a2=AS.atoms[np.logical_and(AS.atoms.resids>=TM2[0],AS.atoms.resids<=TM2[1])]
a6=AS.atoms[np.logical_and(AS.atoms.resids>=TM6[0],AS.atoms.resids<=TM6[1])]

d26a=np.sqrt(((a2.positions.mean(0)-a6.positions.mean(0))**2).sum())

#%% Analyze particular states?



titles=['PCA:1,2/W276:2','PCA:1,2/W276:3','PCA:0/W276:1','PCA:3/W276:1,2','PCA:5/W276:1','PCA:1/W276:2']
t0=np.array([9400,20000,36000,44000,59318,14000])
tf=t0+6000

#Plot where we analyze
fig,ax=plt.subplots(2,1)
ax[0].scatter(np.arange(l[-1])/1000,ECC.PCA.Cluster.state,s=1)
stateW276=np.concatenate([s for s in states.values()])
ax[1].scatter(np.arange(l[-1])/1000,stateW276,s=1)
for a in ax:
    a.set_ylim(a.get_ylim())
    for k,(t,title) in enumerate(zip(t0,titles)):
        color=[c for c in plt.get_cmap('tab10')(k)]
        color[-1]=.25
        a.fill_between([t/1000,t/1000+6],y1=[a.get_ylim()[0],a.get_ylim()[0]],
                       y2=[a.get_ylim()[1],a.get_ylim()[1]],color=color)
    for l0 in l[1:-1]:
        a.plot([l0/1000,l0/1000],a.get_ylim(),color='black')
        
for t,title in zip(t0,titles):
    ax[1].text(t/1000+3,1.5,title.replace('/',' / '),rotation=90,horizontalalignment='center',
               verticalalignment='center')
ax[0].set_ylabel('PCA state')
ax[1].set_ylabel('W276 state')
ax[1].set_xlabel(r't / $\mu$s')
fig.set_size_inches([11.6,5.5])
fig.tight_layout()
fig.savefig(os.path.join(folder,'PCA_chunks.png'))
        
        

proj1=pyDR.Project('Projects/SC_PCA_states',create=True)

# Perform detector / iRED analyses
md_dir='/Volumes/My Book/Y1/apo'

topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k+1}.xtc') for k in range(3)]

select=pyDR.MolSelect(topo,trajs,project=proj1)
select.select_bond('sidechain')
select.traj.step=1

for t00,tf0,title in zip(t0,tf,titles):
    select.traj.t0=t00
    select.traj.tf=tf0
    pyDR.md2data(select,rank=1)
    pyDR.md2iRED(select).iRED2data(rank=1)
    proj1[-2].source.additional_info=title
    proj1[-1].source.additional_info=title

proj1.detect.r_no_opt(10)
proj1.fit()
proj1['no_opt'].detect.r_auto(6)
proj1['no_opt'].fit().opt2dist(rhoz_cleanup=True)
proj1['opt_fit']['iREDmode'].modes2bonds()

proj1.save()

proj2=pyDR.Project('Projects/BB_PCA_states',create=True)
titles=['PCA:1,2/W276:2','PCA:1,2/W276:3','PCA:0/W276:1','PCA:3/W276:1,2','PCA:4,5/W276:1','PCA:1/W276:2']
t0=np.array([9400,20000,36000,44000,59318])
tf=t0+6000


md_dir='/Volumes/My Book/Y1/apo'

topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k+1}.xtc') for k in range(3)]

select=pyDR.MolSelect(topo,trajs,project=proj2)
select.select_bond('N')
select.traj.step=1

for t00,tf0,title in zip(t0,tf,titles):
    select.traj.t0=t00
    select.traj.tf=tf0
    pyDR.md2data(select,rank=1)
    pyDR.md2iRED(select).iRED2data(rank=1)
    proj2[-2].source.additional_info=title
    proj2[-1].source.additional_info=title

proj2.detect.r_no_opt(10)
proj2.fit()
proj2['no_opt'].detect.r_auto(6)
proj2['no_opt'].fit().opt2dist(rhoz_cleanup=True)
proj2['opt_fit']['iREDmode'].modes2bonds()

proj2.save()

#%% Analyze full trajectories
select=pyDR.MolSelect(topo,trajs,project=proj1)
select.select_bond('sidechain')
select.traj.step=1

pyDR.md2data(select,rank=1)
pyDR.md2iRED(select,rank=1).iRED2data()
proj1[-2].source.additional_info='full_traj'
proj1[-1].source.additional_info='full_traj'
proj1['full_traj'].detect.r_no_opt(12)
proj1['full_traj'].fit()
proj1['full_traj']['no_opt'].detect.r_auto(7)
proj1['full_traj']['no_opt'].fit().opt2dist(rhoz_cleanup=True)


select=pyDR.MolSelect(topo,trajs,project=proj2)
select.select_bond('15N')
select.traj.step=1

pyDR.md2data(select,rank=1)
pyDR.md2iRED(select,rank=1).iRED2data()
proj2[-2].source.additional_info='full_traj'
proj2[-1].source.additional_info='full_traj'
proj2['full_traj'].detect.r_no_opt(12)
proj2['full_traj'].fit()
proj2['full_traj']['no_opt'].detect.r_auto(7)
proj2['full_traj']['no_opt'].fit().opt2dist(rhoz_cleanup=True)


#%% Entropy-based CC
select=pyDR.MolSelect(topo,trajs,project=proj1)
select.traj.step=1

commands=['turn x -90','turn y 110','view','lighting full','~ribbon']

for t00,tf0 in zip(t0,tf):
    select.traj.t0=t00
    select.traj.tf=tf0
    select.molsys.make_pdb(ti=0)
    select.select_bond('15N')
    
    ECC=EntropyCC(select)
    proj1.chimera.close()
    ECC.CCchimera(indexCC=ECC.resids==276)
    proj1.chimera.command_line(commands)
    proj1.chimera.savefig(os.path.join(folder,f'ECC_t0_{t00}.png'),
                          options='transparentBackground True',overwrite=True)
    
#%% Total ECC figure
select=pyDR.MolSelect(topo,trajs,project=proj1)
select.traj.step=1
select.select_bond('15N')

ECC=EntropyCC(select)
proj1.chimera.close()
ECC.CCchimera(indexCC=ECC.resids==276)
proj1.chimera.command_line(commands)
proj1.chimera.savefig(os.path.join(folder,f'ECC_total.png'),
                      options='transparentBackground True')

#%% Total iRED CC
proj1.chimera.close()
d=proj1['iREDbond']['full_traj']
d.CCchimera(rho_index=[5],indexCC=d.label==276,scaling=1)
proj1.chimera.command_line(commands)
proj1.chimera.command_line(['ribbon','~show @H*','~sel'])
proj1.chimera.savefig(os.path.join(folder,f'CC_SC_total.png'),options='transparentBackground True',overwrite=True)

proj2.chimera.close()
d=proj2['iREDbond']['full_traj']
d.CCchimera(rho_index=[5],indexCC=d.label==276,scaling=1)
proj2.chimera.command_line(commands)
proj1.chimera.command_line(['~sel'])
proj2.chimera.savefig(os.path.join(folder,f'CC_BB_total.png'),options='transparentBackground True',overwrite=True)