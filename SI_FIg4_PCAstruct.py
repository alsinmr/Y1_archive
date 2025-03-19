#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 13:00:17 2024

@author: albertsmith
"""

from ClusterAna import pca
import os
import matplotlib.pyplot as plt

fig_dir='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/FigureParts/PCAstructure'

proj=pca.project
pca.traj[0]

for k in range(6):
    proj.chimera.close()
    pca.Hist.chimera(n=k,std=2)
    color=','.join([f'{c*100}' for c in plt.get_cmap('tab10')(3)[:-1]])
    proj.chimera.command_line(['~show','ribbon',f'color #2 {color}',
                               'turn y 180','turn z 40','view'])
    proj.chimera.savefig(os.path.join(fig_dir,f'PC{k}_v1.png'),
                         options='transparentBackground True',overwrite=True)
    proj.chimera.command_line(['turn x 90','view'])
    proj.chimera.savefig(os.path.join(fig_dir,f'PC{k}_v2.png'),
                         options='transparentBackground True',overwrite=True)
    