# -*- coding: utf-8 -*-




import pyDR
import os


proj=pyDR.Project('Projects/bound_states',create=True)

step=1

t0=[5900,21200,39700]
length=10800
titles=[f'state{k}' for k in range(len(t0))]


md_dir_as='/Volumes/My Book/Y1/NPY'

topo=os.path.join(md_dir_as,'Y1R_NPY.pdb')
trajs=[os.path.join(md_dir_as,f'run{k+1}.proc.xtc') for k in range(3)]

select=pyDR.MolSelect(topo,trajs,project=proj)
select.traj.step=step
select.select_bond('N')


for k,(t00,title) in enumerate(zip(t0,titles)):
    select.traj.t0=t00
    select.traj.tf=t00+length
    
    pyDR.md2data(select)
    pyDR.md2iRED(select,rank=1).iRED2data()
    
    proj[-2].source.additional_info=title
    proj[-1].source.additional_info=title
    
proj.detect.r_no_opt(10)
proj.fit().detect.r_auto(6)
proj['no_opt'].fit().opt2dist(rhoz_cleanup=True).modes2bonds()

proj.save()