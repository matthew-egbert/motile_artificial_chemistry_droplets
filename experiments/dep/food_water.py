from fr_chemistry import FRChemistry
from molecule import Molecule
from reaction import Reaction
import numpy as np

# oh = Molecule(atom_string='oh', S=0.0, D=0.1, H=0.1)   # pH proxy
# x =  Molecule(atom_string='x',  S=10.0, D=1.0)         # amphiphile
# xx = Molecule(atom_string='xx', S=0.0, D=1.0, H=0.0)   # precursor

a = Molecule(atom_string='a', S=0.0, D=0.5, H=0.0)   # pH proxy
b = Molecule(atom_string='b', S=0.0, D=0.5, H=0.0)   # pH proxy

D = 0.1
y = Molecule(atom_string='y', S=5.0,  D=1.0)
x = Molecule(atom_string='x', S=0.0,  D=D)
w = Molecule(atom_string='w', S=0.0,  D=D)

class ExperimentOSC(object):
    """Two food routes. Oscillates between two food sources. Might as well
be the same food source, AFAICT.

    """
    def __init__(self,model):
        self.model = model

    def reset(self):
        self.model.inds[0].x = 1.
        self.model.inds[0].y = -1.
        self.reset_environment(self.model.env)

    def iterate(self):
        env = self.model.env
        r = 1.0
        amnt = 0.1
        env.modify(-r,0,self.chem.molecules['a'],amnt)
        env.modify(+r,0,self.chem.molecules['b'],amnt)
        # self.model.inds[0].y = -2
        # self.model.inds[0].x = 0.0

    def reset_environment(self,env):
        def g(mesh,xx,yy):
            x,y = mesh
            d = (np.sqrt((xx-x)**2 + (yy-y)**2) < 1.0) * 1.0
            return d

        r = 1.0
        amnt = 0.5*0.5
        a_i = env.get_molecule_index(self.chem.molecules['a'])
        env.c[:,:,a_i] = g(env.coords,-r,0)*amnt
        b_i = env.get_molecule_index(self.chem.molecules['b'])
        env.c[:,:,b_i] = g(env.coords, r,0)*amnt

    def get_chemistry(self):
        class chem(FRChemistry):
            def __init__(self,model):
                self.model = model
                super().__init__()#subsets)

            def reset(self):
                self.subsets = {
                    'env' : [a,b]
                }

                #self.define_molecules([w,x,y,a,b])
                self.define_molecules([w,x,y,b,a])
                self.change_conc_for_all_testtubes(w,1E-1)
                self.change_conc_for_all_testtubes(x,1E-1)
                self.change_conc_for_all_testtubes(a,1E-1)
                self.add_reaction(Reaction([x,a],[x,x],1.,0.0))
                self.add_reaction(Reaction([w,b],[w,w],1.,0.0))
                self.add_reaction(Reaction([x],[y],0.01,0.0))
                self.add_reaction(Reaction([w],[y],0.01,0.0))
                self.add_reaction(Reaction([y],[],0.2,0.0))

                self.add_reaction(Reaction([a,b],[],0.5,0.0))

                
                self.generate_ode_fn()
                self.graph()

            def iterate(self):
                pass


        self.chem = chem(self.model)
        return self.chem
    

# a = Molecule(atom_string='a', S=0.0, D=0.1, H=0.1)   # pH proxy
# b = Molecule(atom_string='b', S=0.00, D=0.1, H=0.1)   # pH proxy
# x =  Molecule(atom_string='x',  S=50.0, D=1.0)         # amphiphile
# xx = Molecule(atom_string='xx', S=0.0, D=1.0, H=0.0)   # precursor


# class ExperimentAB(object):
#     """Basic MOD. Internal fuel that is not replenished. Capable of
#     moving up a gradient.

#     """
#     def __init__(self,model):
#         self.model = model

#     def reset(self):
#         self.model.inds[0].x = -1.
#         self.model.inds[0].y = -1.
#         self.reset_environment(self.model.env)

#     def reset_environment(self,env):
#         def g(mesh,xx,yy):
#             x,y = mesh
#             d = (np.sqrt((xx-x)**2 + (yy-y)**2) < 1.0) * 1.0
#             return d

#         a_i = env.get_molecule_index(self.chem.molecules['a'])
#         env.c[:,:,a_i] = g(env.coords,-1,0)
#         b_i = env.get_molecule_index(self.chem.molecules['b'])
#         env.c[:,:,b_i] = g(env.coords,1,0)

#     def get_chemistry(self):
#         class chem(FRChemistry):
#             def __init__(self,model):
#                 self.model = model
#                 super().__init__()#subsets)

#             def reset(self):
#                 self.subsets = {
#                     'env' : [a,b]
#                 }

#                 self.define_molecules([x,xx,a,b])
#                 self.change_conc_for_all_testtubes(xx,1E-2)
#                 self.add_reaction(Reaction([xx,a],[x,x],0.05,0.0))
#                 self.add_reaction(Reaction([x,b],[xx],0.05,0.0))
#                 # self.add_reaction(Reaction([x],[],0.05,0.0))

#                 self.generate_ode_fn()
#                 self.graph()

#             def iterate(self):
#                 pass


#         self.chem = chem(self.model)
#         return self.chem
