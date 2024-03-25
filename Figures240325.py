#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 11:38:23 2024

@author: albertsmith
"""

import pyDR
import matplotlib.pyplot as plt
import os
import numpy as np

#Projects
HN=pyDR.Project('Projects/basic')
SC=pyDR.Project('Projects/sidechain')

#Data location in projects
titles=['o6:MD:AvOb_state1_0:apo2','o6:MD:AvOb_state2_0:apo1','o6:MD:AvOb_state3_0:apo1']
CCtitles=[title.replace('MD','IREDBOND') for title in titles]

filename=os.path.join('/Users/albertsmith/Documents/Dynamics/Y1_GPCR/Report240325','{}.png')

nd=6

subHN=sum(HN[title] for title in titles) #Subprojects
subSC=sum(SC[title] for title in titles)

#%% Tryptophan amplitudes vs. state




# Setup figure
fig=plt.figure()
ax0=fig.add_subplot(7,1,1)
axHN=[fig.add_subplot(7,2,k) for k in range(3,15,2)]
axSC=[fig.add_subplot(7,2,k) for k in range(4,16,2)]



width=.25
cmap=plt.get_cmap('tab10')
hatches=[None,'xxx','ooo']

# HN motion
i=np.array([sel.resnames[0]=='TRP' for sel in subHN[0].select.sel1],dtype=bool) #Where are tryptophans

for k,(d,hatch) in enumerate(zip(subHN,hatches)):
    for m,(R,a) in enumerate(zip(d.R[i].T,axHN)):
        a.bar(np.arange(i.sum())+(k-1)*width,R,color=cmap(m),width=width,hatch=hatch)
        a.set_xticks(np.arange(i.sum()))
        a.set_xticklabels('')
        a.set_ylim([0,.5])
        a.set_ylabel(fr'$\rho_{m}^{{(\theta,S)}}$',color=cmap(m))

axHN[0].text(0,.35,'Backbone H-N motion')        
axHN[-1].set_ylim([0,1])
fig.set_size_inches([6.75,10.25])
axHN[-1].set_xticklabels(subHN[0].label[i],rotation=90)
axHN[1].legend(('State 1','State 2','State 3'))



# Sidechain motion
i=np.array([sel.resnames[0]=='TRP' for sel in subSC[0].select.sel1],dtype=bool) #Where are tryptophans

for k,(d,hatch) in enumerate(zip(subSC,hatches)):
    for m,(R,a) in enumerate(zip(d.R[i].T,axSC)):
        a.bar(np.arange(i.sum())+(k-1)*width,R,color=cmap(m),width=width,hatch=hatch)
        a.set_xticks(np.arange(i.sum()))
        a.set_xticklabels('')
        a.set_ylim([0,.5])
  
axSC[0].text(0,.35,'Sidechain motion')
axSC[-1].set_ylim([0,1])
axSC[-1].set_xticklabels(subSC[0].label[i],rotation=90)
fig.tight_layout()

        
subHN[0].sens.plot_rhoz(ax=ax0) #Plot sensitivity

fig.savefig(filename.format('trp_motion'))

#%% Just 276W
fig,ax=plt.subplots(2,1)
i=subHN[0].label==276
for m in range(nd):
    for k,(d,hatch) in enumerate(zip(subHN,hatches)):
        ax[0].bar(m+(k-1)*width,d.R[i,m],color=cmap(m),width=width,hatch=hatch)
ax[0].set_xticks(np.arange(nd))
ax[0].set_ylim([0,.12])
ax[0].set_ylabel(r'$\rho_n^{(\theta,S)}$')
ax[0].set_title('HN motion for 276W')

i=subSC[0].label==276
for m in range(nd):
    for k,(d,hatch) in enumerate(zip(subSC,hatches)):
        ax[1].bar(m+(k-1)*width,d.R[i,m],color=cmap(m),width=width,hatch=hatch)
ax[1].set_xticks(np.arange(nd))
ax[1].set_ylim([0,.4])
ax[1].set_ylabel(r'$\rho_n^{(\theta,S)}$')
ax[1].set_title('Sidechain motion for 276W')
ax[1].legend(('State 1','State 2','State 3'))
fig.set_size_inches([5.5,6])
fig.tight_layout()
fig.savefig(filename.format('276W_motion'))

#%% All motion
po=subHN.plot()
po.fig.set_size_inches([7.5,11])
for a in po.ax[:-1]:a.set_ylim([0,.8])
po.ax[1].legend(labels=('State 1','State 2','State 3'),
                handles=[hdl[0][0] for hdl in po.hdls[1]],loc='upper right')
po.ax_sens.set_xlim([-12,-3])
po.show_tc()
subHN.savefig(filename.format('HN_all_motion'))

po=subSC.plot()
po.fig.set_size_inches([7.5,11])
for a in po.ax[:-1]:a.set_ylim([0,.8])
po.ax[1].legend(labels=('State 1','State 2','State 3'),
                handles=[hdl[0][0] for hdl in po.hdls[1]],loc='upper right')
po.ax_sens.set_xlim([-12,-3])
po.show_tc()
subSC.savefig(filename.format('sidechain_all_motion'))

#%% ChimeraX plots of the amplitudes of motion
subHN.chimera()
