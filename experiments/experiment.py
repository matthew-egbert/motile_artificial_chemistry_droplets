import numpy as np
import shutil
import os
import pickle
from collections import defaultdict

from fr_chemistry import FRChemistry
from molecule import Molecule
from reaction import Reaction


class Experiment(object):
    def __init__(self,model,title=''):
        self.model = model
        self.title = title
        
        ## subtrial tracking
        self.trial_start_it = 0        ## last_traj_start_it
        self.trial_start_stop_its = [] ## trajs
        self.n_trials = -1     ## was 21
        self.trial_length = 1000 ## iterations
        self.trial_i  = 0  ## traj_i

        self.RECORDING = False
        
        def constant_factory(value):
            return lambda: value
        self.data = defaultdict(constant_factory(list()))

    def set_to_record(self):
        if self.title == '':
            print("Experiment lacks a title! Don't know where to record data.")
            quit()
            
        self.RECORDING = True            
        self.output_path = f'output/{self.title}/'
        shutil.rmtree(self.output_path,ignore_errors=True)
        os.makedirs(self.output_path)
        
    def write_latex(self):
        if self.RECORDING:
            self.model.surface_chemistry.graph()
            shutil.copyfile('graph.png',self.output_path+'graph.png')
            self.model.surface_chemistry.write_eqs_as_latex(self.output_path+'equations.tex')
            self.model.surface_chemistry.write_ode_latex(self.output_path+'odes.tex')
            self.model.surface_chemistry.write_mol_table_as_latex(self.output_path+'molecules.tex')

    def iterate(self):
        if self.model.it == 0 :
            self.upon_trial_start()        

    def upon_trial_end(self):
        ## record start and stop times of trial        
        self.trial_start_stop_its.append((self.trial_start_it,self.model.it))
        self.trial_start_it = self.model.it

        ## if it was the last trial, wrap things up!
        if self.trial_i == self.n_trials-1:
            self.data['trial_start_stops'] = self.trial_start_stop_its
            self.end_experiment()

        self.trial_i += 1

    def upon_trial_start(self):
        pass
    
a = Molecule(atom_string='a', S=25.0, D=0.0)        
class Testing(Experiment):
    def __init__(self,model) :
        super().__init__('testing',model)
        
    
    def get_chemistry(self):
        class chem(FRChemistry):
            def __init__(self,model):
                self.model = model
                super().__init__()#subsets)

            def reset(self):
                self.subsets = {
                    'env' : [],#[h,p]
                }

                self.define_molecules([a])
                self.set_conc_for_all_testtubes(a,1E-4)
                for o in [24] :
                    for sec_i in [-3,-2,-1,0,1,2,3]:
                        self.model.inds[0].sections[sec_i+o].concentrations[a]=0.1

                self.generate_ode_fn()

            def iterate(self):
                pass

        self.chem = chem(self.model)
        return self.chem

    def reset(self):
        self.model.inds[0].x = -0.
        self.model.inds[0].y = -0.
        #self.reset_environment(self.model.env)

    def iterate(self):
        self.model.inds[0].x = 0.
        self.model.inds[0].y = 0.

        #print('total c: ',sum([s.concentrations[a] for s in self.model.inds[0].sections]))


