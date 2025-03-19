# -*- coding: utf-8 -*-

import pyDR
import os
from pyDR import PCA
from pyDR.PCA import PCA,CombinedPCA

step=1

proj=pyDR.Project()


# Load apo data
md_dir='/Volumes/My Book/Y1/apo'

topo=os.path.join(md_dir,'prot.pdb')
traj=os.path.join(md_dir,'apo1.xtc')

select=pyDR.MolSelect(topo,traj,project=proj)
select.traj.step=step
select.select_bond('15N')

apo=PCA(select)
apo.select_atoms('name N C CA')
apo.align_ref='name CA and resid 28-339'
apo.load()

apo.Movie.xtc_log_swp(nframes=900)

apo.Movie.play_xtc()