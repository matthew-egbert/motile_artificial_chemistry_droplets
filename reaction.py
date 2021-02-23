from typing import Mapping, MutableMapping, Sequence, Iterable, List, Set
from molecule import Molecule
from collections import Counter

class Reaction(object):
    def __init__(self,
                 reactants : List[Molecule],
                 products : List[Molecule],
                 kf : float, kb: float):
        self.rs = reactants
        self.ps = products
        self.kf = kf
        self.kb = kb

    def all_molecules(self):
        return self.rs + self.ps
        
    def __repr__(self):
        r_string = ' + '.join([r.atom_string for r in self.rs])
        p_string = ' + '.join([p.atom_string for p in self.ps])
        return f'{r_string} <--> {p_string}'+ \
               f'(kf={self.kf:.2f}; kb={self.kb:.2f}) \n\t'# + \
               #f'({self.r1.G:.2f} + {self.r2.G:.2f} > {self.p1.G:.2f} + {self.p2.G:.2f}'

    def involves(self,m : Molecule) :
        """ Returns true if r is a reactant or product of this rxn."""
        return (m in self.rs) or (m in self.ps)

    def is_catalyst(self,m : Molecule) :
        return (m in self.rs) and (m in self.ps)

    def loses(self,m : Molecule) :
        """returns the number of the passed molecule lost when reaction
        proceeds (rightwards). Does NOT produce negative values when m
        is gained.

        """
        return max(0,self.rs.count(m) - self.ps.count(m))

    def gains(self,m : Molecule) :
        """returns the number of the passed molecule gained when reaction
        proceeds (rightwards). Does NOT produce negative values when m
        is lost.

        """
        return max(0,self.ps.count(m) - self.rs.count(m))
    
    def on_lhs(self,m : Molecule) :
        """ Returns true if r is a reactant or product of this rxn."""
        return m in self.rs

    def on_rhs(self,m : Molecule) :
        """ Returns true if r is a reactant or product of this rxn."""
        return m in self.ps

    def lhs_count(self,m : Molecule) :
        """ Returns the number of times this molecules is reactant. """
        return self.rs.count(m)

    def rhs_count(self,m : Molecule) :
        """ Returns the number of times this molecule is a product."""
        return self.ps.count(m)
    
    def forward_kinetic_term(self,rxn_i=None):
        if rxn_i is None :
            return f'{self.kf}*{"*".join([r.atom_string for r in self.rs])}'
        else :
            return f'\\phi_{rxn_i}*{"*".join([r.atom_string for r in self.rs])}'

    def backward_kinetic_term(self,rxn_i=None):
        if rxn_i is None :
            return f'{self.kb}*{"*".join([p.atom_string for p in self.ps])}'
        else :
            return f'\\beta_{rxn_i}*{"*".join([p.atom_string for p in self.ps])}'

    def is_meaningless(self):
        reactants = Counter(self.rs)
        products = Counter(self.ps)
        return reactants == products

    def dot_id(self):
        r_string = '+'.join([r.atom_string for r in self.rs])
        p_string = '+'.join([p.atom_string for p in self.ps])
        return f'{r_string}<-->{p_string}'

    def __eq__(self, other):
        reactants = Counter(self.rs)
        other_reactants = Counter(other.rs)
        products = Counter(self.ps)
        other_products = Counter(other.ps)

        return (reactants == other_reactants and products == other_products)

    def latex(self):
        reactants = Counter(self.rs)
        products = Counter(self.ps)

        ## display reactants
        s = '$'
        r = []
        for r_k in reactants.keys() :
            rs = ''
            if reactants[r_k] > 1 :
                rs += str(reactants[r_k])
            rs += r_k.atom_string+' '
            r.append(rs)
        s += '+ '.join(r)
        s += '$ & $'
        if self.kf > 0.0 and self.kb > 0.0 :
            s += r'\leftrightarrow '
        elif self.kf > 0.0 :
            s += r'\rightarrow '
        else :
            r'reaction latex() oddity. '

        ## products
        p = []
        for p_k in products.keys() :
            ss = ''
            if products[p_k] > 1 :
                ss += str(products[p_k])
            ss += p_k.atom_string+' '
            p.append(ss)
        s += '+ '.join(p)

        if len(products) == 0 :
            s+= r'\varnothing '

        s += f'$& {self.kf:0.2} & {self.kb}'
            
        s += r'\\' + '\n'
            
        return(s)
