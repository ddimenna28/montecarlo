"""Top-level package for montecarlo."""

import numpy as np
import math      
import copy as cp
import networkx as nx
import matplotlib.pyplot as plt
import random
import scipy
random.seed(2)

from . bitstring import BitString

from . isingham import IsingHamiltonian

from . montecarlopy import MonteCarlo

'''class BitString:
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
'''

'''class IsingHamiltonian:
    def __init__(self, G):
        """Sets the graph and mus value for the next function 'energy'
        
        Math:
            none
            
        Parameters:
            self: self
                no input
            G: Graph
                input graph initialization
            
        Returns:
            none
        """
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
        """Sets the mus for the next function 'compute_average_values'
        
        Math:
            none
            
        Parameters:
            self: self
                no input
            mus: mus
                input mus initialization
            
        Returns:
            none
            """
        self.mus = mus
    
    def compute_average_values(self, T: float):
        """Computes the average thermodynamic values based on the temperature

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
        """
        
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
'''      
            
'''class MonteCarlo:
    """
    Initialize Bitstring in random configuration
    Create loop that runs n_samples times
        Loop through each site of the Bitstring
            Flip site i
            Check if energy is > or < previous config
                If < keep it
                If >, do prob check
                    If fails, flip site back
    Append energy (e) and magnetization (m)
    """
    def __init__(self, ih):
        """Sets the ih values, including N, for 'run' function
        
            Math:
                none
            
            Parameters:
                self: self
                    no input
                ih: IsingHamiltonian
                    input ising interactions object
            
            Returns:
                none
        """

        N = ih.G.number_of_nodes() # Number of nodes
        self.N = N
        self.ih = ih

        """def calculate_magnetization(self): # Copied from IH above
            spin = 0
            N = self.N # Definitions
            bs = self.bitstring
            for i in range(2**N): # Spin EQ
                bs.set_integer_config(i)
                spin = sum(2 * bs.config[j] - 1 for j in range(N))
            return spin
        """
    
    def run(self, T, n_samples, n_burn):
        """Sets the mus for the next function 'compute_average_values'
        
        Math:
            none
            
        Parameters:
            self: self
                no input
            T: temperature
                input temperature for average values
            n_samples: number of nodes for sample
                input the number of samples to loop through
            n_burn: number of loops burned/removed
                input the number of samples to loop through and remove from the data collection at the beginning of testing
            
        Returns:
            list all_energy (contains old or changed energy values within Ising Hamiltonian objects)
            list all_magnetization (contains old or changed magnetization values within Ising Hamiltonian objects)
        """

        bitstring = BitString(self.N) # Init BitString in Random Config
        for i in range(self.N):
            bitstring.config[i] = np.random.choice([0, 1])
        all_energy = [] # Create empty list for energy
        all_magnetization = [] # Create empty list for magnetization
        current_energy = self.ih.energy(bitstring) # Setting current energy

        for j in range(n_samples): # Loop through BS
            for i in range(self.N): # Loop through each site
                bitstring.flip_site(i) # Flipping site i
                new_energy = self.ih.energy(bitstring) # Setting new energy
                if new_energy < current_energy: # Testing if new energy is less than current energy
                    current_energy = new_energy # Keeping change
                else:
                    boltzman = 1/T
                    prob = np.e**(-boltzman*(new_energy-current_energy)) # Setting Probability
                    if random.random() < prob: # Probability check
                        current_energy = new_energy #Keeping Change
                    else:
                        bitstring.flip_site(i) # FLips site back
            if (n_burn > 0):
                n_burn -= 1
            else:
                all_energy.append(current_energy) # Append to history list
                all_magnetization.append(bitstring.on() - bitstring.off()) # Append to history list
        return all_energy, all_magnetization
'''