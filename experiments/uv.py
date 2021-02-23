from fr_chemistry import FRChemistry
from experiments.experiment import Experiment
from molecule import Molecule
from reaction import Reaction
import numpy as np
import shutil
import os


# ## CHAOS
# FF = 0.02
# kk = 0.047

# DU = 0.010 #* DX2
# DV = 0.005 #* DX2

## DIVIDING
DU = 0.03
DV = 0.005 
FF = 0.055
kk = 0.062

# WIDTH = 2.0*np.pi 
# N_DISCRETIZATIONS = 48
# DX = WIDTH / N_DISCRETIZATIONS
# DX2 = DX**2 

U_DRIVE = FF
U_DEG_K = FF
V_DEG_K = FF+kk

# print(f'U_DRIVE: {U_DRIVE}')
# print(f'U_DEG_K: {U_DEG_K}')
# print(f'V_DEG_K: {V_DEG_K}')

speed_up = 20.0

u = Molecule(atom_string='u', S=0.0, D=DU, H=0.0, X=1.0*speed_up) 
v = Molecule(atom_string='v', S=0.005, D=DV, H=0.0)   

class GrayScottExperiment(Experiment):
    """Basic MOD. Internal fuel that is not replenished. Capable of
    moving up a gradient.

    """
    def __init__(self,model,title='grayscott'):
        super().__init__(model,title)

    def reset(self):
        self.model.inds[0].x = 0.
        self.model.inds[0].y = 0.

        ks = speed_up
        self.model.inds[0].K_R = ks
        self.model.inds[0].K_F = ks
        self.model.inds[0].K_M = ks
        self.reset_environment(self.model.env)

        u_i = self.model.env.get_molecule_index(self.chem.molecules['u'])
        self.model.env.c[:,:,u_i] = U_DRIVE

        def clicked_at(cell_i,cell_j,rx,ry,touch):
            rx = (rx-0.5)*2*self.model.env.petri_r
            ry = (ry-0.5)*2*self.model.env.petri_r
            ## x and y are cell i and j
            r = int(np.sqrt(self.model.env.R / 10))
            if touch.button == 'left':
                d = U_DRIVE*1.2
            else :
                d = 0.0
            self.model.env.c[(self.model.env.coords[0]-rx)**2 +
                             (self.model.env.coords[1]-ry)**2 < r,0] = d

        self.model.env.clicked_at = clicked_at

        
    def iterate(self):
        if self.model.it == 1:
            self.write_latex()
        
        # u_i = self.model.env.get_molecule_index(self.chem.molecules['u'])
        # self.model.env.c[:,:,u_i] = U_DRIVE
        u_i = self.model.env.get_molecule_index(self.chem.molecules['u'])
        self.model.env.c[:,:,u_i] += 0.05*(U_DRIVE - self.model.env.c[:,:,u_i])

    def reset_environment(self,env):
        pass

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
                self.set_conc_for_all_testtubes(v,0.0)
                self.set_conc_for_all_testtubes(u,1.0) 
                self.add_reaction(Reaction([u,v,v],[v,v,v],1.0,0.0)) # 8.0

                self.add_reaction(Reaction([u],[],U_DEG_K,0.0))
                self.add_reaction(Reaction([v],[],V_DEG_K,0.0))

                for o in [0,24] :
                    r = 2
                    for sec_i in range(-r,r+1):
                        self.model.inds[0].sections[sec_i+o].concentrations[v]=1.0
                
                self.generate_ode_fn()
                self.graph()
                
        self.chem = chem(self.model)
        return self.chem
