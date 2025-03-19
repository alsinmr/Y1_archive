#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 12:15:38 2024

@author: albertsmith
"""

import pyDR
import os


#%% Apo W276 state-specific analysis
md_dir='/Volumes/My Book/Y1/apo'
topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k+1}.xtc') for k in range(3)]

proj=pyDR.Project('Projects/apo',create=True)
states=[1,1,1,2,2,2,3]
repeats=[0,1,2,0,1,2,0]
traj_index=[1,2,2,0,0,1,0]
starts=[10000,0,6700,5700,12400,1250,1180,19160]
length=6700

if len(proj['iREDbond'])==0:
    for start,traj,state,repeat in zip(starts,traj_index,states,repeats):
        select=pyDR.MolSelect(topo,trajs[traj],project=proj)
        select.traj.t0=start
        select.traj.tf=start+length
        select.select_bond('N')
        
        pyDR.md2data(select)
        pyDR.md2iRED(select,rank=1).iRED2data()
        
        title=f'W276_state{state}_rep{repeat}'
        proj[-2].source.additional_info=title
        proj[-1].source.additional_info=title
        
    proj.detect.r_no_opt(10)
    proj.fit().detect.r_auto(6)
    proj['no_opt'].fit().opt2dist(rhoz_cleanup=True).modes2bonds()
    
    proj.save()

if len(proj['.+AvOb'])==0:
    from pyDR.misc.Averaging import avgDataObjs
    
    for sub in [proj['opt_fit']['MD'],proj['opt_fit']['iREDbond']]:
        for state in ['state1','state2','state3']:
            avgDataObjs(sub[f'.+{state}'])
            proj[-1].source.additional_info=proj[-1].source.additional_info.rsplit('_',maxsplit=1)[0]
        
proj.save()

#%% NPY cluster-specific analysis
md_dir_as='/Volumes/My Book/Y1/NPY'
topo=os.path.join(md_dir_as,'Y1R_NPY.pdb')
trajs=[os.path.join(md_dir_as,f'run{k+1}.proc.xtc') for k in range(3)]

proj=pyDR.Project('Projects/npy_bound',create=True)
states=[0,1,1,2,3]
repeats=[0,0,1,0,0]
traj_index=[2,1,1,2,0]
starts=[0,3000,9700,10000,10000]
length=6700

if len(proj['iREDbond'])==0:
    for start,traj,state,repeat in zip(starts,traj_index,states,repeats):
        select=pyDR.MolSelect(topo,trajs[traj],project=proj)
        select.traj.t0=start
        select.traj.tf=start+length
        select.select_bond('N')
        
        pyDR.md2data(select)
        pyDR.md2iRED(select,rank=1).iRED2data()
        
        title=f'cluster{state}_rep{repeat}'
        proj[-2].source.additional_info=title
        proj[-1].source.additional_info=title
        
    proj.detect.r_no_opt(10)
    proj.fit().detect.r_auto(6)
    proj['no_opt'].fit().opt2dist(rhoz_cleanup=True).modes2bonds()
    
    proj.save()

if len(proj['.+AvOb'])==0:
    from pyDR.misc.Averaging import avgDataObjs
    
    for sub in [proj['opt_fit']['MD'],proj['opt_fit']['iREDbond']]:
        for state in [f'cluster{k}' for k in range(4)]:
            avgDataObjs(sub[f'.+{state}'])
            proj[-1].source.additional_info=proj[-1].source.additional_info.rsplit('_',maxsplit=1)[0]
        
proj.save()


#%% Apo full length analylisis
proj=pyDR.Project('Projects/full_length_ana',create=True)

length=10000
starts=[6359,6307,5492,16359,16307]
tis=[0,1,2,0,1]

md_dir='/Volumes/My Book/Y1/apo'

topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k+1}.xtc') for k in range(3)]

for k,(start,ti) in enumerate(zip(starts,tis)):
    select=pyDR.MolSelect(topo,trajs[ti],project=proj)
    select.select_bond('N')
    
    select.traj.t0=start
    select.traj.tf=start+length
    
    pyDR.md2data(select)
    pyDR.md2iRED(select,rank=1).iRED2data()
    
    title=f'apo_rep{k}'
    proj[-2].source.additional_info=title
    proj[-1].source.additional_info=title


#%% NPY analysis
tis=[0,1,2]
starts=[6736,7063,6942]



md_dir_as='/Volumes/My Book/Y1/NPY'

topo=os.path.join(md_dir_as,'Y1R_NPY.pdb')
trajs=[os.path.join(md_dir_as,f'run{k+1}.proc.xtc') for k in range(3)]

for k,(start,ti) in enumerate(zip(starts,tis)):
    select=pyDR.MolSelect(topo,trajs[ti],project=proj)
    select.select_bond('N')
    
    select.traj.t0=start
    select.traj.tf=start+length
    
    pyDR.md2data(select)
    pyDR.md2iRED(select,rank=1).iRED2data()
    
    title=f'NPY_rep{k}'
    proj[-2].source.additional_info=title
    proj[-1].source.additional_info=title
    


#%% Load antagonist-bound data
tis=[0,1,2,0,1,2]
starts=[2047,1876,1862,12047,11876,11862]

md_dir_as='/Volumes/My Book/Y1/Antagonist'

topo=os.path.join(md_dir_as,'Y1R_UR_MK299.pdb')
trajs=[os.path.join(md_dir_as,f'run{k+1}.xtc') for k in range(3)]

for k,(start,ti) in enumerate(zip(starts,tis)):
    select=pyDR.MolSelect(topo,trajs[ti],project=proj)
    select.select_bond('N')
    
    select.traj.t0=start
    select.traj.tf=start+length
    
    pyDR.md2data(select)
    pyDR.md2iRED(select,rank=1).iRED2data()
    
    title=f'anta_rep{k}'
    proj[-2].source.additional_info=title
    proj[-1].source.additional_info=title


#%% Load NPY/Gi-bound data
md_dir_as='/Volumes/My Book/Y1/NPY_Gi'

tis=[0,1,2]
starts=[501,701,1000]

topo=os.path.join(md_dir_as,'Y1R_NPY_Gi.pdb')
trajs=[os.path.join(md_dir_as,f'run{k+1}.proc.xtc') for k in range(3)]

for k,(start,ti) in enumerate(zip(starts,tis)):
    select=pyDR.MolSelect(topo,trajs[ti],project=proj)
    select.select_bond('N')
    
    select.traj.t0=start
    select.traj.tf=start+length
    
    pyDR.md2data(select)
    pyDR.md2iRED(select,rank=1).iRED2data()
    
    title=f'Gi_rep{k}'
    proj[-2].source.additional_info=title
    proj[-1].source.additional_info=title
    

#%% Fit all data from NYP, antagonist, Gi boud

proj.detect.r_no_opt(10)
proj.save()

proj.fit().detect.r_auto(6)
proj['no_opt'].fit().opt2dist(rhoz_cleanup=True).modes2bonds()

proj.save()

#%% Average repetitions
    
if len(proj['.+AvOb'])==0:
    from pyDR.misc.Averaging import avgDataObjs
    
    for sub in [proj['opt_fit']['MD'],proj['opt_fit']['iREDbond']]:
        for traj in ['apo','NPY','anta','Gi']:
            avgDataObjs(sub[f'.+{traj}'])
            proj[-1].source.additional_info=proj[-1].source.additional_info.rsplit('_',maxsplit=1)[0]
            
proj.save()