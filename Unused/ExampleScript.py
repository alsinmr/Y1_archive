# -*- coding: utf-8 -*-

import pyDR
import os
from pyDR.PCA.PCAclean import PCA
from CombinedPCA import CombinePCA

step=1

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


#Basic stuff
apo.Lambda  #PCA eigenvalues
apo.PC      #List of principle components
apo.PCxyz   #Reshaping of above to separate x,y and z components
apo.PCamp   #Time dependence of the principle components
apo.ref_pos #Reference position used for alignment



# This plots the histogram
apo.Hist.plot()

# This connects to ChimeraX and shows structures
apo.Hist.hist2struct(cmap='Reds',from_traj=False)  #from_traj takes a specific frame