# -*- coding: utf-8 -*-



import pyDR
import os
from pyDR.PCA import PCA,CombinedPCA

step=10

TRP_res=[106,148,163,276,288]

proj=pyDR.Project()


# Load apo data
md_dir='/Volumes/My Book/Y1/apo'

topo=os.path.join(md_dir,'prot.pdb')
traj=os.path.join(md_dir,f'apo1.xtc')

select=pyDR.MolSelect(topo,traj,project=proj)
select.traj.step=step
select.select_bond('N')

pca=PCA(select)
pca.select_atoms('protein and not name H*')
pca.align_ref='name CA and resid 28-339'
pca.load()


pca.Movie.xtc_log_swp(nframes=900)
pca.Movie.play_xtc()