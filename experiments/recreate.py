from pylab import *
from experiments.experiment import Experiment
from fr_chemistry import FRChemistry
from molecule import Molecule
from reaction import Reaction
import pickle
from sty import fg, bg, ef, rs

factor = 0.5#0.13089969*0.25

x_factor = 5.0


h = Molecule(atom_string='H', S=0.0,   D=0.1*factor, H=0.1, X = 1.0*x_factor)   # pH proxy
a = Molecule(atom_string='A', S=50.0,  D=1.0*factor)                   # amphiphile
p = Molecule(atom_string='P', S=0.0,   D=1.0*factor, H=0.0, X = 0.25*x_factor)  # precursor


def g(xx,yy,mesh,sigma=2.5):
    x,y = mesh
    d = np.sqrt((xx-x)**2 + (yy-y)**2)
    a = 1.0   # max
    c = sigma # stddev
    g = a*np.exp(-(d)**2/(2.*c**2))
    return g

h_mag = 0.0001*0.2
h_x   = -1.5
p_mag = 0.0001*0.2
p_x   = 1.5

## MOD initial position
init_pos = (0.0,-4.0)

class RecreateExperiment(Experiment):

    def __init__(self,model,
                 refuelling=False,
                 dynamic_env=False,
                 motion_allowed=False,
                 title_extension=''):

        title = '__'.join([
            {True : 'refuelling',
             False: 'no_refuelling',}[refuelling],
            {True : 'dynamic_env',
             False: 'fixed_env',}[dynamic_env],
            {True : 'motion',
             False: 'no_motion'}[motion_allowed]])

        title+=title_extension
        
        self.refuelling = refuelling
        self.dynamic_env = dynamic_env
        self.motion_allowed = motion_allowed

        ## no refuelling exp (1A-1D)
        if (self.refuelling == False and
            self.dynamic_env == False) :
            self.final_it = 800
        elif (self.refuelling == False and
            self.dynamic_env == True) :            
            self.final_it = 800                

        ## refuelling exp (2A-2D)
        elif (self.refuelling == True and
              self.dynamic_env == False) :
            self.final_it = 1000
        elif (self.refuelling == True and
              self.dynamic_env == True) :
            self.final_it = 5000
        
        super().__init__(model,title)

    def reset(self) :
        self.model.inds[0].x = init_pos[0]
        self.model.inds[0].y = init_pos[1]
        self.model.inds[0].K_R = 1.0
        self.model.inds[0].K_F = 1.0
        self.model.inds[0].K_M = 1.0
        self.reset_environment()

    def reset_environment(self):
        env = self.model.env
        h_i = env.get_molecule_index(self.chem.molecules['H'])
        p_i = env.get_molecule_index(self.chem.molecules['P'])

        env.c[:,:,h_i] = h_mag*g(h_x,0.0,env.coords)
        if self.refuelling :
            env.c[:,:,p_i] = p_mag*g(p_x,0.0,env.coords)
        
        env.mu[h_i] = 0.005
        env.mu[p_i] = 0.005

    def iterate(self):
        super().iterate()
        if self.model.it == 1 :
            self.write_latex()

        if self.motion_allowed == False :
            self.model.inds[0].x = init_pos[0]
            self.model.inds[0].y = init_pos[1]

        if self.dynamic_env == False:
            self.reset_environment()
            
        ## check if time to end entire experiment
        if (self.model.it) == self.final_it:
            print(f'apparently {self.model.it} = {self.final_it}')
            self.end_experiment()

    def end_experiment(self):
        if self.RECORDING:
            self.model.plot_scene(show_time=True)
            plt.savefig(self.output_path+'/scene.png',dpi=300)

            self.data['ind'] = self.model.inds[0].get_data()
            self.data['ts']  = self.model.ts
            self.data['env'] = {}
            for env_chem in self.model.surface_chemistry.subsets['env']:
                ec = env_chem
                ec_i = self.model.env.get_molecule_index(ec)
                self.data['env'][ec.atom_string] = np.array(self.model.env.c[:,:,ec_i])
            self.data['env']['mask'] = self.model.env.mask

            with open(self.output_path+'/data.p', 'wb') as handle:
                pickle.dump(dict(self.data), handle, protocol=pickle.HIGHEST_PROTOCOL)
            quit()

    def get_chemistry(self):
        class chem(FRChemistry):
            def __init__(self,model):
                self.model = model
                super().__init__()#subsets)

            def reset(self):
                self.subsets = {
                    'env' : [h,p]
                }

                self.define_molecules([a,p,h])
                self.set_conc_for_all_testtubes(a,0.9E-4)
                self.set_conc_for_all_testtubes(p,0.9E-4)
                self.set_conc_for_all_testtubes(h,0.0)
                self.add_reaction(Reaction([p,h],[a,a],10000.0,0.0))

                self.generate_ode_fn()

        self.chem = chem(self.model)
        return self.chem
        
class ExperimentVaryChi(RecreateExperiment):#ExperimentRefuellingFixed):
    """I repeat the static environment refuelling experiment, each time
    changing the exchange rate of p, to show that the behaviour is
    responding to the needs of the agent.

    """
    def __init__(self,model,title_extension='__vary_chi'):
        super().__init__(model,refuelling=True,
                         dynamic_env=False,
                         motion_allowed=True,
                         title_extension=title_extension)
        
        self.trial_length = 1500
        #self.chis = np.geomspace(0.05,5.0,self.n_trials)#[0.05,0.1,0.25,0.5,1.0]#
        self.chis = np.geomspace(0.1,12.5,11)#[0.05,0.1,0.25,0.5,1.0]#
        self.n_trials = len(self.chis)
        self.data['chis'] = self.chis

        ## diables typical end of trial mechanism, allowing subtrial
        ## system to decide when to finish
        self.final_it = -1 

    def iterate(self):
        super().iterate()
        
        if self.model.inds[0].y > 0.0 :
            self.upon_trial_end()
            self.upon_trial_start()

    def upon_trial_start(self) :
        ### reset trial
        chi = self.chis[self.trial_i]
        p.X = chi
        self.model.inds[0].x = init_pos[0]
        self.model.inds[0].y = init_pos[1]
        self.model.surface_chemistry.reset()
        print(bg.da_blue+f'({self.model.it}) starting trial {self.trial_i} | Setting: chi:{p.X}'+bg.rs)

    def upon_trial_end(self) :
        super().upon_trial_end()
        print(bg.da_blue+f'({self.model.it}) trial {self.trial_i} recorded: {self.trial_start_stop_its[-1]} | values were: chi:{p.X}  final_x: {self.model.inds[0].x}'+bg.rs)






