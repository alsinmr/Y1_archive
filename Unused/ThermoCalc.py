# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

class ThermoCalc():
    k=1.380649e-23
    R=8.31446261815324
    def __init__(self,pca,nbins=10,T_K=310,states=None):
        """
        Calculates Thermodynamic parameters a list of states in a principle
        component analysis based on binning of the PCA. Returned values are
        relative energies/entropies, not absolute values.
        
        PCs will be used only up until the bin width becomes larger than
        3x the standard deviation of the  PC.
        
        Caution should be taken with this approach. The free energy is easily
        biased by undersampling, and by the starting state of the simulation.

        Parameters
        ----------
        pca : PCA object
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        self.pca=pca
        self.nbins=nbins
        self.states=states
        self.T_K=T_K
        
        self._bins=None
        self._PCstate=None
        
    @property
    def nbins(self):
        return self._nbins
    @nbins.setter
    def nbins(self,nbins):
        self._bins=None
        self._nbins=nbins
        self._PCstate=None
        
    @property
    def states(self):
        return self._states
    @states.setter
    def states(self,states):
        states=np.array(states)
        assert np.issubdtype(states.dtype,np.integer),"States must be integers"
        assert states.size==len(self.pca.traj),"states must be the same length as the trajectory"
        self._PCstate=None
        self._states=states
        
    @property
    def n_states(self):
        return len(np.unique(self.states))
        
    @property
    def bins(self) -> np.array:
        """
        Bins to use for PCA binning

        Returns
        -------
        np.array

        """
        if self._bins is None:
            PCmax=np.abs(self.pca.PCamp).max()
            self._bins=np.linspace(-PCmax,PCmax,self.nbins+1)
        return self._bins
    
    @property
    def bin_width(self) -> float:
        """
        Width of a bin

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return self.bins[1]-self.bins[0]
    
    @property
    def nPCs(self) -> int:
        """
        Number of principal components to bin for the entropy estimate

        Returns
        -------
        int

        """
        return np.sum(np.sqrt(self.pca.Lambda)*2>self.bin_width)
    
    
    @property
    def PCstate(self):
        if self._PCstate is None:
            state0=np.sum([self.pca.PCamp[:self.nPCs]>bin0 for bin0 in self.bins[1:]],axis=0).astype('uint64')
            self._PCstate=state0
            # state=np.zeros(state0.shape[-1],dtype=np.uint64)
            # for k in range(self.nPCs):
            #     state+=state0[k]*(self.nbins**k)
            # self._PCstate=state
                
        return self._PCstate
        
    @property
    def S(self):
        S=[]
        for state0 in np.arange(self.states.max()+1):
            count=np.unique(self.PCstate[:,state0==self.states],return_counts=True,axis=-1)[1]
            p=count/(state0==self.states).sum()
            S.append(-self.R*(p*np.log(p)).sum())
        return np.array(S)
    
    @property
    def population(self):
        p=[]
        for state0 in np.arange(self.states.max()+1):
            p.append((state0==self.states).sum()/len(self.states))
        return np.array(p)
    
    @property
    def PCpopulation(self):
        count=np.unique(self.PCstate,return_counts=True,axis=-1)[1]
        p=count/len(self.pca.traj)
        return p
    
    @property
    def DelG(self):
        p=self.population
        G=-np.log(p)*self.T_K*self.R
        G-=G.min()
        return G
    
    @property
    def DelH(self):
        out=self.DelG+self.T_K*self.S
        out-=out.min()
        return out
    
    def plot(self):
        fig,ax=plt.subplots(4,1,figsize=[4.8,6.6])
        cmap=plt.get_cmap('tab10')
        ax[0].bar(np.arange(self.n_states),self.population,color=cmap(0))
        ax[0].set_ylabel('Population')
        ax[1].bar(np.arange(self.n_states),self.DelG/1e3,color=cmap(1))
        ax[1].set_ylabel(r'$\Delta$ G / kJ*mol$^{-1}}$')
        ax[2].bar(np.arange(self.n_states),self.S,color=cmap(2))
        ax[2].set_ylabel(r'S / J*mol$^{-1}$K$^{-1}$')
        ax[3].bar(np.arange(self.n_states),self.DelH/1e3,color=cmap(3))
        ax[3].set_ylabel(r'$\Delta$ H / kJ*mol$^{-1}}$')
        fig.tight_layout()
        
        for a in ax:
            a.set_xticks(np.arange(self.n_states))
            a.set_xlabel('State')
        
        return ax
        
            
        
            
            
        
        
        
    
        
        