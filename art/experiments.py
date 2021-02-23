from fr_chemistry import FRChemistry
from molecule import Molecule
from reaction import Reaction
import numpy as np

class Experiment1(object):
    def __init__(self,model):
        self.model = model

    def reset_environment(self,env):
        env.c[:,:,0] = (1.0-np.linspace(0,1,env.R))*2.0
        env.c[:,:,1] = 0.0*np.linspace(0,1,env.R)

    def get_chemistry(self):
        class chem(FRChemistry):
            def __init__(self,model):
                self.model = model
                super().__init__()#subsets)

            def reset(self):
                oh = Molecule(atom_string='oh', S=0.0,  D=0.)   # pH proxy
                x =  Molecule(atom_string='x',  S=0.1, D=0.01) # amphiphile
                xx = Molecule(atom_string='xx', S=0.0,  D=0.99)  # precursor        

                self.subsets = {
                    'env' : [oh]#,xx]
                }

                self.define_molecules([x,xx,oh])
                self.change_conc_for_all_testtubes(xx,1E0)
                self.add_reaction(Reaction([xx,oh],[x,x],0.05,0.0))
                self.add_reaction(Reaction([x],[],0.1,0.0))

                self.generate_ode_fn()
                self.graph()

            def iterate(self):
                pass

        return chem(self.model)


class Experiment2(object):
    def __init__(self,model):
        self.model = model

    def reset_environment(self,env):
        env.c[:,:,0] = (1.0-np.linspace(0,1,env.R))*1.0
        env.c[:,:,1] = 0.5

    def get_chemistry(self):
        class chem(FRChemistry):
            def __init__(self,model):
                self.model = model
                super().__init__()#subsets)

            def reset(self):
                oh = Molecule(atom_string='oh', S=0.0,  D=0.)   # pH proxy
                x =  Molecule(atom_string='x',  S=1.0, D=0.01) # amphiphile
                xx = Molecule(atom_string='xx', S=0.0,  D=0.99)  # precursor        

                self.subsets = {
                    'env' : [oh,xx]
                }

                self.define_molecules([x,xx,oh])
                self.change_conc_for_all_testtubes(xx,1E0)
                self.add_reaction(Reaction([xx,oh],[x,x],0.05,0.0))
                self.add_reaction(Reaction([x],[],0.1,0.0))

                self.generate_ode_fn()
                self.graph()

            def iterate(self):
                pass

        return chem(self.model)
    

