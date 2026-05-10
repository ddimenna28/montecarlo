import numpy as np
import math      
import copy as cp
import networkx as nx
import matplotlib.pyplot as plt
import random
import scipy
random.seed(2)

from . bitstring import BitString

class IsingHamiltonian:
    def __init__(self, G):
        '''Sets the graph and mus value for the next function 'energy'
        
        Math:
            none
            
        Parameters:
            self: self
                no input
            G: Graph
                input graph initialization
            
        Returns:
            none
        '''
        self.G = G
        self.mus = np.zeros(G.number_of_nodes())

    def energy(self, bs: BitString):
        """Compute energy of configuration, `bs`

            .. math::
                E = \\left<\\hat{H}\\right>

        Parameters
        ----------
        bs   : Bitstring
            input configuration
        G    : Graph
            input graph defining the Hamiltonian
        Returns
        -------
        energy  : float
            Energy of the input configuration
        """
        N = len(bs)
        G = self.G
        E = 0.0
        for e in G.edges:
            i,j = e
            weight = G.edges[e]['weight']
            si = 2 * bs.config[i] - 1
            sj = 2 * bs.config[j] - 1
            E += weight * si * sj
        for i in range(N):
            E += self.mus[i] * (2 * bs.config[i] - 1)
        return E
    
    def set_mu(self, mus):
        '''Sets the mus for the next function 'compute_average_values'
        
        Math:
            none
            
        Parameters:
            self: self
                no input
            mus: mus
                input mus initialization
            
        Returns:
            none'''
        self.mus = mus
    
    def compute_average_values(self, T: float):
        '''Computes the average thermodynamic values based on the temperature

        Math:
            <E> = Sum over index a of (E(a)*Probability(a))
            <M> = Sum over index a of (M(a)*Probability(a))
            HC = (<E^2>-<E>^2)*T^(-2)
            MS = (<M^2>-<M>^2)*T^(-1)
            Contains averages
        
        Parameters:
            self: self
                no input
            T: temperature (float)
                input temperature that all graph energies are dependent on
        
        Returns:
            Average Energy
            Magnetization
            Variance of E
            Variance of M
        '''
        
        G = self.G
        bs = BitString(G.number_of_nodes())
        N = len(bs)

        E  = 0.0 #Average Energy
        M  = 0.0 #Magnetization
        Z  = 0.0 
        EE = 0.0 #Variance of E
        MM = 0.0 #Variance of M
        MS = 0.0
        
        boltzman = 1/T #boltzman constant

        # Finding Z
        for i in range(2**N):
            bs.set_integer_config(i) # Finding E(alpha)
            Z += np.e**(-boltzman*self.energy(bs))
        
        # Finding Average Energy
        for i in range(2**N):
            bs.set_integer_config(i)
            E += self.energy(bs)*((np.e**(-boltzman*self.energy(bs)))/(Z))
        
        # Finding Heat Capacity
        for i in range(2**N):
            bs.set_integer_config(i)
            EE += ((self.energy(bs))**2)*((np.e**(-boltzman*self.energy(bs)))/(Z))
        HC = (EE-E**2) / (T**2)

        # Finding Magnetization
        for i in range(2**N):
            bs.set_integer_config(i)
            spin = sum(2 * bs.config[j] - 1 for j in range(N))
            M += spin*((np.e**(-boltzman*self.energy(bs)))/(Z))

        # Finding Magnetization Susceptibility
        for i in range(2**N):
            bs.set_integer_config(i)
            MM += ((sum(2 * bs.config[j] - 1 for j in range(N)))**2)*((np.e**(-boltzman*self.energy(bs)))/(Z))
        MS = (MM - M**2) / T
        
        
        
        return E, M, HC, MS