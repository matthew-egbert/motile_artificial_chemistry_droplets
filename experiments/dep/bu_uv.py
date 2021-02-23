from fr_chemistry import FRChemistry
from experiments.experiment import Experiment
from molecule import Molecule
from reaction import Reaction
import numpy as np
import shutil
import os

MAG = 5.0
FF = 0.02 * MAG
kk = 1.0*0.047 * MAG
U_DRIVE = FF
U_DEG_K = FF
V_DEG_K = FF+kk


print(U_DRIVE)
# quit()

# oh = Molecule(atom_string='oh', S=0.0, D=0.1, H=0.1)   # pH proxy
# x =  Molecule(atom_string='x',  S=10.0, D=1.0)         # amphiphile
# xx = Molecule(atom_string='xx', S=0.0, D=1.0, H=0.0)   # precursor

v = Molecule(atom_string='v', S=0.5, D=0.0001, H=0.0)   # pH proxy
#f = Molecule(atom_string='f',  S=0.0, D=0.5, H=0.0)   # pH proxy

u = Molecule(atom_string='u', S=0.0,  D=0.0001, H=0.5,X=10.0)
# y = Molecule(atom_string='y', S=5.0,  D=1.0)
# w = Molecule(atom_string='w', S=0.0,  D=D)

# SETUP RVIT TO BE ABLE TO CHANGE CHEMICAL PARAMETERS DYNAMICALLY

class ExperimentWaste(Experiment):
    """Basic MOD. Internal fuel that is not replenished. Capable of
    moving up a gradient.

    """
    def __init__(self,model,title='uv'):
        super().__init__('testing',model)

    def reset(self):
        self.model.inds[0].x = -0.
        self.model.inds[0].y = -4.
        self.reset_environment(self.model.env)

        def clicked_at(cell_i,cell_j,rx,ry,touch):
            rx = (rx-0.5)*2*self.model.env.petri_r
            ry = (ry-0.5)*2*self.model.env.petri_r
            ## x and y are cell i and j
            r = int(np.sqrt(self.model.env.R / 10))
            if touch.button == 'left':
                d = U_DRIVE*5
            else :
                d = 0.0
            self.model.env.c[(self.model.env.coords[0]-rx)**2 +
                             (self.model.env.coords[1]-ry)**2 < r,0] = d

        self.model.env.clicked_at = clicked_at


    def iterate(self):
        if self.model.it == 1:
            self.write_latex()

        env = self.model.env
        r = 1.0
        amnt = 0.1
        #env.modify(-r,0,self.chem.molecules['u'],-amnt)
        u_i = env.get_molecule_index(self.chem.molecules['u'])
        # env.c[:,:,u_i] += self.model.p1#g(env.coords,-r,0)*amnt
        
        env.c[:,:,u_i] = U_DRIVE*self.model.p1


    def reset_environment(self,env):
        def g(mesh,xx,yy):
            x,y = mesh
            d = (np.sqrt((xx-x)**2 + (yy-y)**2) < 1.0) * 1.0
            return d

        u_i = env.get_molecule_index(self.chem.molecules['u'])
        env.c[:,:,u_i] = U_DRIVE#self.model.p1 * 100.0# self.model.p1 # 0.4
        # env.c[:int(env.R/2),:,u_i] = 0.7# self.model.p1 # 0.4
        env.mu[u_i] = 0.1

    def get_chemistry(self):
        class chem(FRChemistry):
            def __init__(self,model):
                self.model = model
                super().__init__()#subsets)

            def reset(self):
                self.subsets = {
                    'env' : [u]
                }

                
                self.define_molecules([u,v])
                self.change_conc_for_all_testtubes(v,0.04)
                self.change_conc_for_all_testtubes(u,U_DRIVE) 
                self.add_reaction(Reaction([u,v,v],[v,v,v],LAMBDA,0.0)) # 8.0

                k = 1.0
                self.add_reaction(Reaction([u],[],U_DEG_K,0.0))
                self.add_reaction(Reaction([v],[],V_DEG_K,0.0))

                for o in range(0,48,16) :
                    for sec_i in [-3,-2,-1,0,1,2,3]:
                        self.model.inds[0].sections[sec_i+o].concentrations[v]=0.0
                        self.model.inds[0].sections[sec_i+o].concentrations[u]=U_DRIVE
                # for i in range(24) :
                #     self.model.inds[0].sections[i].concentrations[v]=0.0
                
                self.generate_ode_fn()
                self.graph()
                
            def iterate(self):
                pass


        self.chem = chem(self.model)
        return self.chem
