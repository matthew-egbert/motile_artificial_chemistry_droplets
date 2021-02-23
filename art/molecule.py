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
                 #H : float = None, ## hydrophobicity
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
            print(self.atom_string)
            self.atom_count = atom_count

        self.length = len(self.atom_string)
        ## G is free energy
        if G != None:
            self.G = G
        else :
            low,high = G_range
            self.G = np.random.rand()*(high-low)+low
            print(' %f < %f < %f ' %(low,self.G,high))

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
            
        ## Take any remaining kwargs and store them as attributes.
        ## This allows special molecules (e.g. surfactants) to have
        ## extra values associated with them, such as surface-tension
        for key, value in kwargs.items():
            setattr(self, key, value)

    def dot_id(self):
        s = self.atom_string
        s += '\n(%f)' %(self.G)
        return s

    def __repr__(self):
        return f'{self.atom_string} : G={self.G:.2f} S={self.S:.2f}'
        

    def __eq__(self, other):        
        return other.atom_string==self.atom_string

    def __hash__(self):
        return self.atom_string.__hash__()
