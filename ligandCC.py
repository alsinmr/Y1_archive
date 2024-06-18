# -*- coding: utf-8 -*-

import pyDR
import os


md_dir='/Volumes/My Book/Y1/NPY'

topo=os.path.join(md_dir,'Y1R_NPY.pdb')
trajs=[os.path.join(md_dir,f'run{k}.proc.xtc') for k in range(1,4)]

proj=pyDR.Project('Projects/ligandCC',create=True)

select=pyDR.MolSelect(topo,trajs,project=proj)
select.select_bond('N')

pyDR.md2iRED(select,rank=1).iRED2data()

proj1=pyDR.Project('Projects/BB_PCA_states')['iREDbond']