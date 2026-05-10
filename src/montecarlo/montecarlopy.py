import numpy as np
import math      
import copy as cp
import networkx as nx
import matplotlib.pyplot as plt
import random
import scipy
random.seed(2)

class MonteCarlo:
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