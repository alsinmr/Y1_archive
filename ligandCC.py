# -*- coding: utf-8 -*-

import pyDR
import os

proj=pyDR.Project('Projects/ligandCC',create=True)



#%% Backbone
md_dir='/Volumes/My Book/Y1/NPY'

topo=os.path.join(md_dir,'Y1R_NPY.pdb')
trajs=[os.path.join(md_dir,f'run{k}.proc.xtc') for k in range(1,4)]



select=pyDR.MolSelect(topo,trajs,project=proj)
select.select_bond('N')

pyDR.md2iRED(select,rank=1).iRED2data()

proj.detect.r_no_opt(12)
proj.fit().detect.r_auto(7)
proj['no_opt'].fit().opt2dist(rhoz_cleanup=True).modes2bonds()


md_dir='/Volumes/My Book/Y1/NPY_Gi'
topo=os.path.join(md_dir,'Y1R_NPY_Gi.pdb')
trajs=[os.path.join(md_dir,f'run{k}.proc.xtc') for k in range(1,4)]

select=pyDR.MolSelect(topo,trajs,project=proj)
select.select_bond('N')

pyDR.md2iRED(select,rank=1).iRED2data()

proj[-1].detect.r_no_opt(12)
proj[-1].fit().detect.r_auto(7)
proj[-1].fit().opt2dist(rhoz_cleanup=True).modes2bonds()

proj.save()

#%% Sidechain
md_dir='/Volumes/My Book/Y1/NPY'

topo=os.path.join(md_dir,'Y1R_NPY.pdb')
trajs=[os.path.join(md_dir,f'run{k}.proc.xtc') for k in range(1,4)]

select=pyDR.MolSelect(topo,trajs,project=proj)
select.select_bond('sidechain')

pyDR.md2iRED(select,rank=1).iRED2data()
proj[-1].source.additional_info='SC_npy'

proj[-1].detect.r_no_opt(12)
proj[-1].fit().detect.r_auto(7)
proj[-1].fit().opt2dist(rhoz_cleanup=True).modes2bonds()


md_dir='/Volumes/My Book/Y1/NPY_Gi'
topo=os.path.join(md_dir,'Y1R_NPY_Gi.pdb')
trajs=[os.path.join(md_dir,f'run{k}.proc.xtc') for k in range(1,4)]

select=pyDR.MolSelect(topo,trajs,project=proj)
select.select_bond('sidechain')

pyDR.md2iRED(select,rank=1).iRED2data()
proj[-1].source.additional_info='SC_npy_Gi'

proj[-1].detect.r_no_opt(12)
proj[-1].fit().detect.r_auto(7)
proj[-1].fit().opt2dist(rhoz_cleanup=True).modes2bonds()

proj.save()

#%% PCA analysis