#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 14:49:16 2025

@author: albertsmith
"""

import os
from ClusterAna import ECC,pca
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pyDR


proj=pca.project

# Figure directory
fig_dir='/Users/albertsmith/Documents/DynamicsLaptop/Y1_GPCR/Figures/Figure2'

color=[plt.get_cmap('turbo').resampled(5)(k) for k in range(1,5)]
cmapTrp=ListedColormap(color)



#%% Figure: PCA Backbone CC


#W276 rotameric clustering (Ramachandran plot, part C)
index=np.argmax(ECC.resids==276)
ax=ECC.plotChi(index,cmap=cmapTrp)[0]
fig=ax.figure
fig.set_size_inches([5,5])
ax.set_aspect('equal')
fig.savefig(os.path.join(fig_dir,'W276_Ramachandran.png'))


# Plot backbone/sidechain cross-correlation (ECC, part A)
ax=ECC.plotCCpca()
fig=ax.figure
ax.set_xticks(np.arange(0,len(ECC.resi),10))
fig.set_size_inches([12.5,4])
fig.tight_layout()

fig.savefig(os.path.join(fig_dir,'sidechain_pca_CC.pdf'))


# 3D Structures (part B)
ECC.CCchimera(indexCC='PCA')
ECC.project.chimera.savefig(os.path.join(fig_dir,'sidechain_pca_CC_3D.png'),options='transparentBackground True')
