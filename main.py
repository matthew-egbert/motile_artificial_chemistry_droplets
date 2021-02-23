from pylab import *

from kivy.config import Config
Config.set('kivy','log_level','info')
# Config.set('kivy','log_enable', '0')
Config.write()
from kivy.clock import Clock

from rvit.core import init_rvit,rvit_reconnect

from individual import Individual
from molecule import Molecule
from reaction import Reaction
from environment import Environment
from event_logger import EventLogger

from experiments.experiment import Testing
from experiments.tricol import *
from experiments.recreate import *
from experiments.uv import *
# from experiments.testing import *
# from experiments.food_water import *

#FernandoRoweChemistry, ByHandChemistry
from sty import fg, bg, ef, rs

import sys,time

class Drawing(object):
    def __init__(self,model):
        self.model = model
        self.ind_xs = np.array([i.x for i in self.model.inds])
        self.ind_ys = np.array([i.y for i in self.model.inds])

        self.N_SECTIONS = self.model.inds[0].N_SECTIONS

        ## FLOW LINES
        self.f_lines = np.zeros((self.model.N_INDS,
                                 self.N_SECTIONS*2,
                                 2))

        ## CONCENTRATION RINGS
        self.MAX_CONCENTRATION_RINGS = 10
        self.c_rings = np.zeros((self.model.N_INDS,
                                 self.MAX_CONCENTRATION_RINGS,
                                 self.N_SECTIONS,
                                 2))

        from itertools import chain
        wrap_indices = np.array( list(chain(range(1,self.N_SECTIONS*2),[0])) )
        n_verts = self.MAX_CONCENTRATION_RINGS*self.model.N_INDS*self.N_SECTIONS
        self.c_ring_indices = [repeat(arange(self.N_SECTIONS)+n,2)[wrap_indices]
                               for n in range(0,n_verts,self.N_SECTIONS)]
        self.c_ring_indices = ravel(self.c_ring_indices)
        self.c_ring_colors  = np.zeros((self.model.N_INDS,
                                        self.MAX_CONCENTRATION_RINGS,
                                        self.N_SECTIONS))

        for ind_i,ind in enumerate(self.model.inds):
            for mol_i in range(self.MAX_CONCENTRATION_RINGS):
                for section_i in range(self.N_SECTIONS):
                    self.c_ring_colors[ind_i,mol_i,section_i] = \
                        float(mol_i) / self.MAX_CONCENTRATION_RINGS

    def iterate(self):
        self.ind_xs[:] = [i.x for i in self.model.inds]
        self.ind_ys[:] = [i.y for i in self.model.inds]

        ## concentration rings
        mols_to_draw = self.model.surface_chemistry.molecules.values()
        mols_to_draw = list(mols_to_draw)[:self.MAX_CONCENTRATION_RINGS]
        for ind_i,ind in enumerate(self.model.inds):
            for mol_i,mol in enumerate(mols_to_draw):
                for section_i in range(self.N_SECTIONS):
                    conc = np.log10(ind.sections[section_i].concentrations[mol])
                    clip_min = -4
                    clip_max = 1
                    conc = max(min(clip_max,conc),clip_min)
                    conc = 2.0*(conc-clip_min)/(clip_max-clip_min)
                    x = ind.x + cos(ind.thetas[section_i])*conc
                    y = ind.y + sin(ind.thetas[section_i])*conc
                    self.c_rings[ind_i,mol_i,section_i,:] = x,y

        for ind_i,ind in enumerate(self.model.inds):
            for section_i in range(self.N_SECTIONS):
                theta = ind.thetas[section_i]+np.pi/self.N_SECTIONS
                r = ind.R
                x = ind.x + cos(theta)*r
                y = ind.y + sin(theta)*r
                self.f_lines[ind_i,section_i*2,:] = x,y
                flow_amnt = 10000.0*ind.convective_flux[section_i]
                dx = cos(theta+np.pi/2)*flow_amnt
                dy = sin(theta+np.pi/2)*flow_amnt
                self.f_lines[ind_i,section_i*2+1,:] = x+dx,y+dy


class Model(object):
    def __init__(self,run_experiment=None):

        self.PETRI_R = 10.0
        self.DT = 0.05
        self.ts = []

        ## UI parameters
        self.p1 = 1.0
        self.p2 = 0.2

        seed = np.random.randint(100)
        #seed = 61
        print(f'SEED is {seed}')
        np.random.seed(seed)

        if run_experiment == None :
            self.set_experiment(0)
        else:
            print(bg.blue + f'=== Running Experiment {run_experiment} ===' + bg.rs)
            self.set_experiment(run_experiment,recording=True)

        def iterate(arg):
            self.iterate()
            self.drawing.iterate()
            self.ts.append(self.t)
            self.t += self.DT
            if self.it % 500 == 0 :
                print(f'time: {self.t:.2f} it: {self.it}')

        Clock.schedule_interval(iterate, 1.0 / 100.0)
        self.drawing = Drawing(self)

    def set_experiment(self,exp_i,recording=False):
        experiments = [
            ## 0: Exp 1A
            RecreateExperiment(self,
                               refuelling=False,
                               dynamic_env=False,
                               motion_allowed=True),
            ## 1: Exp 1B
            RecreateExperiment(self,
                               refuelling=False,
                               dynamic_env=False,
                               motion_allowed=False),
            ## 2: Exp 1C
            RecreateExperiment(self,
                               refuelling=False,
                               dynamic_env=True,
                               motion_allowed=False),
            ## 3: Exp 1D
            RecreateExperiment(self,
                               refuelling=False,
                               dynamic_env=True,
                               motion_allowed=True),

            ## 4: Exp 2A
            RecreateExperiment(self,
                               refuelling=True,
                               dynamic_env=False,
                               motion_allowed=True),
            ## 5: Exp 2B
            RecreateExperiment(self,
                               refuelling=True,
                               dynamic_env=False,
                               motion_allowed=False),
            ## 6: Exp 2C
            RecreateExperiment(self,
                               refuelling=True,
                               dynamic_env=True,
                               motion_allowed=False),
            ## 7: Exp 2D
            RecreateExperiment(self,
                               refuelling=True,
                               dynamic_env=True,
                               motion_allowed=True),            

            ## 8
            ExperimentVaryChi(self),

            ## 9
            GrayScottExperiment(self),
        ]

        self.experiment = experiments[exp_i]
        if recording:            
            self.experiment.set_to_record()
            
        self.t = 0.0
        self.ts = []
        self.it = 0
        self.is_paused = False
        self.surface_chemistry = self.experiment.get_chemistry()

        self.N_INDS = 1
        self.env = Environment(self)
        self.inds = [Individual(self) for _ in range(self.N_INDS)]
        self.surface_chemistry.reset()
        self.env.reset()
        self.experiment.reset()

        rvit_reconnect()

    def pause(self):
        self.is_paused = True

    def play(self):
        self.is_paused = False

    def iterate(self):
        if not self.is_paused:
            self.env.iterate()
            self.experiment.iterate()
            for i in self.inds:
                i.iterate()
            self.it += 1

    def trigger_cascade(self):
        n_rxn = len(self.surface_chemistry.reactions)
        while( len(self.surface_chemistry.reactions) == n_rxn):
            print('------CASCADE TRIGGERED MANUALLY----')
            self.surface_chemistry.generate_avalanche()
        self.surface_chemistry.graph()

    def plot_scene(self,show_plot=False,**kwargs):
        fig = figure(figsize=(8,8))
        ax = gca()
        r = self.env.petri_r
        ind = self.inds[0]
        c = self.env.c[:,:,0]
        c[self.env.mask] = c.max()
        imshow(c,cmap='gray',extent=[-r,r,-r,r],origin='lower',interpolation='bilinear')

        petri = plt.Circle((0., 0.), r,
                            color='None',ec='k', lw=10, clip_on=True)
        ax.add_artist(petri)

        mod = plt.Circle((ind.x, ind.y), ind.R,
                            color='r', clip_on=False)
        ax.add_artist(mod)
        ax.axis('off')
        tight_layout()
        k = 1.05*r
        xlim(-k,k)
        ylim(-k,k)
        plot(ind.x_h[10:],ind.y_h[10:],'r--')

        if hasattr(kwargs,'show_time') :
            title(f'{self.t:0.2}')

        # ## plot polar concentrations
        # ax = plt.subplot(gs[0,0], projection='polar')
        # cs = np.array([[section.concentrations[m]
        #                 for m in self.model.surface_chemistry.molecules.values()]
        #                   for section in self.sections])
        # ts = hstack([self.thetas,self.thetas[0]])
        # ts = vstack([ts,]*shape(cs)[1])
        # cs = vstack([cs[:,:],cs[0,:][newaxis]])
        # ax.plot(ts.T, log(cs))
        # # ax.set_rmax(2)
        # # ax.set_rticks([0.5, 1, 1.5, 2])  # less radial ticks
        # ax.set_rlabel_position(-22.5)  # get radial labels away from plotted line
        # #ax.grid(True)

        if show_plot:
            show()


if __name__ == '__main__':
    if len(sys.argv)>1 :
        m = Model(run_experiment=int(sys.argv[1]))
        r = 1
        init_rvit(m,rvit_file='rvit.kv',window_size=(1200/r,800/r))
    else :
        m = Model()
        r = 1
        init_rvit(m,rvit_file='rvit.kv',window_size=(1320,880))
        # init_rvit(m,rvit_file='rvit.kv',window_size=(1200/r,800/r))
    # m.event_logger.close()
