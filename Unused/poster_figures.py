#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 13:40:13 2025

@author: albertsmith
"""

index=ECC.resids==276
i=(ECC.state[index]==0).flatten()



fig,ax=plt.subplots()
ax.scatter(pca.pcas[-1].PCamp[0],pca.pcas[-1].PCamp[1],c=plt.get_cmap('tab10')(1),s=.1)
ax.scatter(pca.pcas[-2].PCamp[0],pca.pcas[-2].PCamp[1],c=plt.get_cmap('tab10')(2),s=.1)
ax.scatter(pca.pcas[0].PCamp[0][i],pca.pcas[0].PCamp[1][i],c=plt.get_cmap('tab10')(0),s=.1)
ax.set_aspect(1)
ax.set_xlim([-150,150])
ax.set_ylim([-150,150])
fig.savefig('/Volumes/No Name/Y1_Figures/Poster1.png',transparent=True)


i=(ECC.state[index]==1).flatten()
fig,ax=plt.subplots()
ax.scatter(pca.pcas[-1].PCamp[0],pca.pcas[-1].PCamp[1],c=plt.get_cmap('tab10')(1),s=.1)
ax.scatter(pca.pcas[-2].PCamp[0],pca.pcas[-2].PCamp[1],c=plt.get_cmap('tab10')(2),s=.1)
ax.scatter(pca.pcas[0].PCamp[0][i],pca.pcas[0].PCamp[1][i],c=plt.get_cmap('tab20c')(1),s=.1)
ax.set_aspect(1)
ax.set_xlim([-150,150])
ax.set_ylim([-150,150])
fig.savefig('/Volumes/No Name/Y1_Figures/Poster2.png',transparent=True)

i=(ECC.state[index]==2).flatten()
fig,ax=plt.subplots()
ax.scatter(pca.pcas[-1].PCamp[0],pca.pcas[-1].PCamp[1],c=plt.get_cmap('tab10')(1),s=.1)
ax.scatter(pca.pcas[-2].PCamp[0],pca.pcas[-2].PCamp[1],c=plt.get_cmap('tab10')(2),s=.1)
ax.scatter(pca.pcas[0].PCamp[0][i],pca.pcas[0].PCamp[1][i],c=plt.get_cmap('tab20c')(3),s=.1)
ax.set_aspect(1)
ax.set_xlim([-150,150])
ax.set_ylim([-150,150])
fig.savefig('/Volumes/No Name/Y1_Figures/Poster3.png',transparent=True)