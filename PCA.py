#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 09:04:26 2024

@author: albertsmith
"""

import pyDR
import os
from pyDR.PCA import PCA,CombinedPCA

step=1

TRP_res=[106,148,163,276,288]

proj=pyDR.Project()


# Load apo data
md_dir='/Volumes/My Book/Y1/apo'

topo=os.path.join(md_dir,'prot.pdb')
trajs=[os.path.join(md_dir,f'apo{k+1}.xtc') for k in range(3)]

select=pyDR.MolSelect(topo,trajs,project=proj)
select.traj.step=step

apo=PCA(select)
apo.select_atoms('name N C CA and resid 28-339')
apo.align_ref='name CA and resid 28-339'
apo.load()




#%% Load NPY-bound data
md_dir_as='/Volumes/My Book/Y1/NPY'

topo=os.path.join(md_dir_as,'Y1R_NPY.pdb')
trajs=[os.path.join(md_dir_as,f'run{k+1}.proc.xtc') for k in range(3)]

select=pyDR.MolSelect(topo,trajs,project=proj)
select.traj.step=step

bound=PCA(select)
bound.select_atoms('name N C CA and segid B')

# Non-standard usage: we copy the reference positions from the apo trajectory
# and use them here so that we can align the two trajectories together.
bound._ref_pos = apo.ref_pos
bound.align_ref='name CA and segid B'
bound.load()

#%% Load antagonist-bound data
md_dir_as='/Volumes/My Book/Y1/Antagonist'

topo=os.path.join(md_dir_as,'Y1R_UR_MK299.pdb')
trajs=[os.path.join(md_dir_as,f'run{k+1}.xtc') for k in range(3)]

select=pyDR.MolSelect(topo,trajs,project=proj)
select.traj.step=step

anta=PCA(select)
anta.select_atoms('name N C CA and segid R and resid 28-339')

# Non-standard usage: we copy the reference positions from the apo trajectory
# and use them here so that we can align the two trajectories together.
anta._ref_pos = apo.ref_pos
anta.align_ref='name CA and segid R and resid 28-339'
anta.load()


#%% Load NPY/Gi-bound data
md_dir_as='/Volumes/My Book/Y1/NPY_Gi'

topo=os.path.join(md_dir_as,'Y1R_NPY_Gi.pdb')
trajs=[os.path.join(md_dir_as,f'run{k+1}.proc.xtc') for k in range(3)]

l_gi=[]
for traj in trajs:
    select=pyDR.MolSelect(topo,traj)
    select.traj.step=step
    l_gi.append(len(select.traj))

select=pyDR.MolSelect(topo,trajs,project=proj)
select.traj.step=step

gi=PCA(select)
gi.select_atoms('name N C CA and segid D and resid 28-339')

# Non-standard usage: we copy the reference positions from the apo trajectory
# and use them here so that we can align the two trajectories together.
gi._ref_pos = apo.ref_pos
gi.align_ref='name CA and segid D and resid 28-339'
gi.load()




All=CombinedPCA(apo,bound,anta,gi)
three=CombinedPCA(apo,bound,anta)
