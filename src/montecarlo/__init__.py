"""Top-level package for montecarlo."""

import numpy as np
import math      
import copy as cp
import networkx as nx
import matplotlib.pyplot as plt
import random
import scipy
random.seed(2)


class BitString:
    """
    Simple class to implement a config of bits
    """
    def __init__(self, N):
        self.N = N
        self.config = np.zeros(N, dtype=int) 

    def __repr__(self):
        out = ""
        for i in self.config:
            out += str(i)
        return out

    def __eq__(self, other):        
        return all(self.config == other.config)
    
    def __len__(self):
        return len(self.config)

    def on(self):
        """
        Return number of bits that are on
        """
        on_bits = 0
        for i in self.config:
            if i == 1:
                on_bits += 1
        return on_bits

    def off(self):
        """
        Return number of bits that are off
        """
        off_bits = self.N
        for i in self.config:
            if i == 1:
                off_bits -= 1
        return off_bits

    def flip_site(self,i):
        """
        Flip the bit at site i
        """
        if self.config[i] == 0:
            self.config[i] = 1
        else:
            self.config[i] = 0
    
    def integer(self):
        """
        Return the decimal integer corresponding to BitString
        """
        decimal = 0
        for i in range(self.N):
            decimal += self.config[i] * pow(2, self.N-1-i)
        return decimal
 

    def set_config(self, s:list[int]):
        """
        Set the config from a list of integers
        """
        for i in range(self.N):
            self.config[i] = s[i]

    def set_integer_config(self, dec:int):
        """
        convert a decimal integer to binary
    
        Parameters
        ----------
        dec    : int
            input integer
            
        Returns
        -------
        Bitconfig
        """
        for i in range(self.N):
            self.config[self.N-1-i] = dec % 2
            dec = dec // 2
        return BitString(self.N)

class IsingHamiltonian:
    def __init__(self, G):
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
        self.mus = mus
    
    def compute_average_values(self, T: float):
        
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
            


