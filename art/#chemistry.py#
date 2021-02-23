from typing import Mapping, MutableMapping, Sequence, Iterable, List, Set
import numpy as np
from matplotlib import pyplot as plt
from collections import defaultdict
import itertools
from graphviz import Digraph
from jinja2 import Template
from scipy.integrate import odeint

from molecule import Molecule
from reaction import Reaction

class ConcentrationHistory(object):
    def __init__(self,mol):
        self.mol = mol
        self.ts = []
        self.cs = []

    def plot(self,no_labels=False,**kwargs):
        if no_labels:
            plt.plot(self.ts,self.cs,**kwargs)
        else :
            plt.plot(self.ts,self.cs,label=self.mol.atom_string,**kwargs)

class TestTube(object):
    """This object describes a location where chemistry takes place. It is
    associated with a single instance of a Chemistry (but many
    TestTubes can be associated with a single Chemistry). It tracks
    concentration histories, and uses the associated Chemistry to
    update the concentrations. Other things that change the
    concentrations (e.g. diffusion between TestTubes) can be applied
    by directly changing the self.concentrations. """

    def __init__(self, chemistry):
        self.chemistry = chemistry
        self.chemistry.associate_testtube(self)

        ## key is a Molecule
        ## value is float indicating concentration
        def constant_factory(value):
            return lambda: value
        self.concentrations = defaultdict(constant_factory(0))

        ## tracks the history of any molecules that appear
        self.conc_h = {}

    def sample(self,t):
        for i,m in enumerate(self.chemistry.molecules.values()):
            if m not in self.conc_h.keys() :
                self.conc_h[m] = ConcentrationHistory(m)
            self.conc_h[m].ts.append(t)
            self.conc_h[m].cs.append(self.concentrations[m])


    def get_concentration_histories(self,include=None,exclude=None,**kwargs):
        histories = {}
        if include == None:
            include = self.chemistry.molecules.values()
        if exclude == None:
            exclude = []
        for mol in self.conc_h.keys() :
            if mol in include and mol not in exclude:
                histories[mol] = self.conc_h[mol]
        return histories

    def plot_concentration_histories(self,include=None,exclude=None,no_labels=False,**kwargs):
        histories = self.get_concentration_histories(include,exclude,**kwargs)

        for mol in histories.keys():
            ## for consistent color plotting, get the index of this molecule
            ## from the chemistry
            color = self.chemistry.get_color_of_molecule(mol)
            self.conc_h[mol].plot(no_labels=no_labels,color=color,**kwargs)


class Chemistry(object):
    """This object describes a set of molecules and a set of rules that
    can be applied to change those molecules. It does not describe any
    concentrations of molecules (for that, see TestTube).

    """

    def get_color_of_molecule(self,mol):
        mol_index = list(self.molecules.keys()).index(mol.atom_string)
        mol_index = mol_index % 10
        color = 'C%d'%(mol_index)
        #print('%s : %s' %(mol,color))
        return color

    def __init__(self, subsets={}):
        """subsets is a dictionary of string -> list[molecules].  These
        subsets allow certain molecules to be treated differently
        (e.g. certain molecules have no intrinsic dynamics, but change
        only due to external ("environmental") reasons. Other sets are
        facilitate certain calculations -- e.g. to sum the concentration
        of certain molecules.

        """

        ## key is a string showing the molecule, e.g.: aab
        ## value is a Molecule
        self.molecules = {}

        self.testtubes = []

        self.subsets = subsets
        for subset in self.subsets.values():
            for m in subset:
                self.define_molecule(m)

        self.reactions = []
        self.testtubes = []

    def associate_testtube(self,testtube:TestTube):
        self.testtubes.append(testtube)

    def global_concentration(self,m:Molecule):
        return sum([tt.concentrations[m] for tt in self.testtubes])            

    def define_molecules(self,mols:List[Molecule]):
        for mol in mols:
            self.define_molecule(mol)
    
    def define_molecule(self,molecule : Molecule):#, concentration=0.0):
        #print('defining molecule %s' %(molecule))
        # if concentration > 0.0:
        #     self.change_conc_for_all_testtubes(molecule,concentration)
        self.molecules[molecule.atom_string] = molecule

    def change_conc_for_all_testtubes(self,molecule : Molecule, concentration):
        for tt in self.testtubes:
            tt.concentrations[molecule] += concentration        
        
    def add_reaction(self,rxn : Reaction):
        if rxn.is_meaningless():
            print('Refusing to add meaningless reaction.')
            return
        self.reactions.append(rxn)

    def graph(self):
        #with open("chemistry.txt", "w") as text_file:
        print(f"{self}", file=open('khemistry.txt','w'))

        ## below this concentration, molecules are not drawn
        c_threshold = 0#1E-15

        dot = Digraph(comment='chem',engine='neato') # neato
        for m in self.molecules.values():
            v = self.global_concentration(m)
            scaled_v = min(0.98,max(0.3,0.15*-np.log10(1E-15+v)))
            color=f'0.5 0.0 {scaled_v}'
            label = f'<<b>{m.atom_string}</b><br/>' +\
                    f'G:{m.G:.2f}<br/>' +\
                    f'[{np.log10(1E-15+v):.2f}]>'
            atts = {
                'fontsize' : '12',
                'fontname' : 'roboto',
            }
            if m in self.subsets['env']:
                dot.node(m.dot_id(),shape='rect',fillcolor='gold',style='filled',
                         label=label,**atts)
            else:
                shape = 'circle'
                if v < c_threshold :
                    shape = 'point'
                    label = ''
                dot.node(m.dot_id(),shape=shape,
                         style='filled',color=color,
                         label=label, **atts)

        for rxn in self.reactions:
            if any([self.global_concentration(m) >=
                    c_threshold for m in rxn.all_molecules()]) :
                dot.node(rxn.dot_id(),shape='point')
                catalysts_added = []
                for m in rxn.all_molecules() :
                    if rxn.is_catalyst(m):
                        if m not in catalysts_added:
                            dot.edge(m.dot_id(),rxn.dot_id(),arrowhead='icurve')
                            catalysts_added.append(m)
                    if rxn.loses(m):
                        ## color="black:black:"
                        dot.edge(rxn.dot_id(),m.dot_id(),arrowhead='inv')
                    if rxn.gains(m):
                        dot.edge(rxn.dot_id(),m.dot_id())

        # print(dot.source)
        dot.format = 'png'
        dot.render('graph')

    def kinetics_rhs(self,m):
        """generates the right-hand-side of an ODE that describes m's
        kinetics.

        """

        s = '0.0 '

        if m in self.subsets['env']:
            return s + ' # env variable'
        
        for r in self.reactions:
            if self.a_constituent_in_this_reaction_exists(r) :
                #print(r)
                for _ in range(r.lhs_count(m)):
                    s+= '-'+r.forward_kinetic_term()
                    if r.kb != 0.0:
                        s+= '+'+r.backward_kinetic_term()
                for _ in range(r.rhs_count(m)):
                    s+= '+'+r.forward_kinetic_term()
                    if r.kb != 0.0:
                        s+= '-'+r.backward_kinetic_term()

        # if m not in self.subsets['env'] :
        #     s += f'-0.1*{m.atom_string}'
        import sympy
        simplified_expr = sympy.sympify(s)
        return simplified_expr

    def a_constituent_in_this_reaction_exists(self,rxn) :
        ## possible optimization here, but could be complicated
        ## need to track which mols actually exist in simulation
        ## and cleverly regenerate the ODEs when new ones appear
        return True

        # this is hacky don't use it
        # for tt in self.testtubes:
        #     for mol in [rxn.r1,rxn.r2,rxn.p1,rxn.p2] :
        #         if tt.concentrations[mol] > 1E-5:
        #             return True
        # return False


    def generate_ode_fn(self):
        """Creates a method `self.ode_fn(y,t)` which calculates the
        instantaneous change in all of the dynamic molecules (i.e. all
        that do not belong to the subset 'env'). This method is then
        added to this instance.

        """

        template = """
def dydt(y,t):
    # load state variables
{{molecules}} = y

    # calculate derivatives
{{derivs}}
    # return derivatives
    return {{d_molecules}}
        """
        t = Template(template)

        all_mols = self.molecules.values()
        dynamic_mols = [m for m in self.molecules.values()]# if m not in self.subsets['env']]
        molecules = '    '+','.join([m.atom_string for m in all_mols])
        d_molecules = ','.join(['d_'+m.atom_string for m in all_mols])
        derivs = '\n'.join('    d_%s = %s'%(m.atom_string,self.kinetics_rhs(m))
                           for m in dynamic_mols)

        fn = t.render(molecules=molecules,
                      derivs=derivs,
                      d_molecules=d_molecules)
        print(fn)
        ns = {}
        exec(fn,ns)
        self.ode_fn = ns['dydt']

    # def simulate(self, testtube, t, ode_fn=None):
    #     """ Simulates the chemistry from some IC. Used to test the generated ode fn."""

    #     y0 = [testtube.concentrations[m] for m in self.molecules.values()]
    #     if ode_fn == None:
    #         ode_fn = self.ode_fn

    #     sol = odeint(ode_fn, y0, t)
    #     n = np.shape(sol)[1]

    #     for i,m in enumerate(self.molecules.values()):
    #         if m not in testtube.conc_h.keys() :
    #             testtube.conc_h[m] = ConcentrationHistory(m)
    #         testtube.conc_h[m].ts.extend(t)
    #         testtube.conc_h[m].cs.extend(sol[:,i])

    #     ## update concentrations
    #     for i,m in enumerate(self.molecules.values()):
    #         testtube.concentrations[m] = sol[-1,i]

    #     return sol,


    def dcdt(self, testtube, dt, ode_fn=None):
        """Calculates how the internal chemistry dynamics will change the
        concentrations. Returns the CHANGE IN concentrations.
        """

        y0 = [testtube.concentrations[m] for m in self.molecules.values()]
        if ode_fn == None:
            ode_fn = self.ode_fn

        sol = odeint(ode_fn, y0, np.linspace(0,dt,10))
        # n = np.shape(sol)[1]

        ## update concentrations
        for i,m in enumerate(self.molecules.values()):
            testtube.concentrations[m] = sol[-1,i]

        assert(np.all(sol[0,:] == y0))
        change = sol[-1,:] - sol[0,:]

        return change

    def plot_concentration_history(self,include=None,exclude=None):
        if include == None:
            include = self.molecules.values()
        if exclude == None:
            exclude = []
        for mol in self.conc_h.keys() :
            if mol in include and mol not in exclude:
                self.conc_h[mol].plot()

    def cull_unconnected_molecules(self):
        """Removes all molecules that do not take part in any reactions AND
        are not in any subsets.

        """
        mols_to_remove = []
        for m in list(self.molecules.values()):
            m_is_unconnected = True
            for r in self.reactions:
                if r.involves(m):
                    m_is_unconnected = False
                    break
            if m_is_unconnected:
                mols_to_remove.append(m)

        for m in mols_to_remove:
            m_is_in_a_subset = False
            for subset in self.subsets.values():
                if m in subset:
                    m_is_in_a_subset = True
            if not m_is_in_a_subset:
                self.molecules.pop(m.atom_string)

    def cull_irrelevant_reactions(self):
        """Experimental. Any reaction that has low reactant and product
        concentration has a chance of being removed."""

        ## define a function f, that maps log10(rxn_cs) to chances of being removed
        ## it is a linear piecemeal function, below `a` it is guaranteed to be
        ## removed and above `b` is is guaranteed to remain. Inbetween values
        ## are interpreted linearly.
        
        # remove_below
        a = -12
        #keep_above
        b = -5

        m = -1.0 / (b-a)
        k = 1.0 - (m*a)

        def f(c):
            if c <= a:
                return 1.0
            if c >= b:
                return 0.0
            return m*c+k

        ## the concentration of all of the molecules
        ## relevant to each reaction in self.reactions
        for rxn in self.reactions:
            rxn_c = np.log10(sum([self.global_concentration(m)
                               for m in rxn.all_molecules()]))
            p_remove = f(rxn_c)
            print(f'{rxn_c} removal chance: {p_remove}')
            if np.random.rand() < p_remove:
                print('REMOVED')
                self.reactions.remove(rxn)
        
        #print(rxn_c)
        # score = 1./(1E-15+c_of_all_mols)
        # score = score/sum(score)
        # print(score)
        #quit()
        # print(f'%{p}  for total conc of: {total_c}')
        # TODO DECIDE ODDS OF REMOVING REACTION
        # if np.random.rand() < p  :
        #     print(f'Removing {rxn}')
        #     self.reactions.remove(rxn)
        print(f'#m: {len(self.molecules)} \t #r:{len(self.reactions)}')
        

    def __repr__(self):
        s = '\n\n'
        s += 'MOLECULES\n'
        s += '=========\n'
        total_volume = 0.0
        for key in self.molecules.keys():
            mol = self.molecules[key]
            #conc = self.concentrations[mol]
            s+= '%s: (%f)' %(key.rjust(10),mol.G)
            #s+= ('%0.4f mM\n'%(conc)).rjust(13)
            #total_volume += conc
        #s+= 'Total volume: %f\n' %(total_volume)

        s += '\n\n'
        s += 'REACTIONS\n'
        s += '=========\n'
        for rxn in self.reactions:
            s+=str(rxn) + '\n'
        return s
