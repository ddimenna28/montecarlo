import numpy as np

class BitString:
    """
    Simple class to implement a config of bits
    """
    def __init__(self, N):
        """
        Initializes configuration
        """
        self.N = N
        self.config = np.zeros(N, dtype=int) 

    def __repr__(self):
        """
        Reports configuration
        """
        out = ""
        for i in self.config:
            out += str(i)
        return out

    def __eq__(self, other):
        """
        Tests configuration equality
        """  
        return all(self.config == other.config)
    
    def __len__(self):
        """
        Returns length of configuration
        """
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