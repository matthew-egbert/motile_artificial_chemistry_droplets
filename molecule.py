import numpy as np
from utils import logger

atoms = ('a','b','c','d')

class Molecule(object):
    def __init__(self,
                 atom_string : str = None,
                 atom_count : list = [],
                 G : float = None,
                 G_range : tuple = (0,1),
                 D : float = None, ## diffusion rate
                 S : float = None, ## surfactancy constant
                 H : float = None, ## hydrophilia
                 X : float = None, ## env/interface exchange rate
                 **kwargs):
        """If atom_count is specified, it indicates the number of each atom
        that should be in the molecule. [2,3,0,0] means 2A, 3B and no
        C or D atoms.

        G :: free energy -- can be specified directly as a float or
        alternatively G_range specifies the minimum and maximum
        allowed values for a flat distribution that it will be
        selected from. Units are kJ/mol.

        """
        global atoms

        if atom_string is not None :
            self.atom_string = atom_string
            self.atom_count = [self.atom_string.count(c) for c in self.atom_string]
        if atom_count != []:
            atoms = [item for a,n in zip(atoms,atom_count) for item in [a,]*n]
            np.random.shuffle(atoms)
            self.atom_string = ''.join(atoms)
            self.atom_count = atom_count

        self.length = len(self.atom_string)
        ## G is free energy
        if G != None:
            self.G = G
        else :
            low,high = G_range
            self.G = np.random.rand()*(high-low)+low

        ## D is diffusion rate constant
        if D != None:
            self.D = D
        else :
            low,high = 0.0,1.0
            self.D = np.random.rand()*(high-low)+low

        ## S is surfactancy constant
        if S != None:
            self.S = S
        else :
            low,high = -0.1,0.1
            self.S = np.random.rand()*(high-low)+low


        ## H is hydrophilia (only relevant for env mols)
        #### 1.0 means will only leave droplet
        #### 0.0 means will only enter droplet
        #### middle values are also allowed
        if H != None:
            self.H = H
        else :
            low,high = 0.0,1.0
            self.H = 0.5


        ## X is surface/env exchange rate
        #### 0.0 means no exchange
        #### middle values are also allowed
        if X != None:
            self.X = X
        else :
            self.X = 1.0
            
        ## Take any remaining kwargs and store them as attributes.
        ## This allows special molecules (e.g. surfactants) to have
        ## extra values associated with them, such as surface-tension
        for key, value in kwargs.items():
            setattr(self, key, value)

    def dot_id(self):
        s = self.atom_string
        s += '\n(%f)' %(self.G)
        return s

    def table_latex(self):
        return r'$'+f'{self.atom_string}'+'$' +f'& {self.D} & {self.S} & {self.H} & {self.X}\\\\' + '\n'

    def __repr__(self):
        return f'{self.atom_string} : G={self.G:.2f} D={self.D:.2f} S={self.S:.2f} H={self.H:.2f}'
    

    def __eq__(self, other):        
        return other.atom_string==self.atom_string

    def __hash__(self):
        return self.atom_string.__hash__()
